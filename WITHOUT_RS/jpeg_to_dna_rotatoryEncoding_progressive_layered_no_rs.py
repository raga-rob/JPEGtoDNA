"""
jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py

Encodes a JPEG file into fixed-length DNA fragments WITHOUT Reed-Solomon error correction:
- Priority packing of earlier JPEG bytes into earlier fragments (leveraging progressive JPEG).
- Optional whitening (XOR with keystream) on fragment body to decorrelate structure.
- Rotating 2-bit->base mapping with per-fragment phase to reduce repeats.
- Keystream-derived, non-periodic padding to reach exact -l bases per fragment.
- Option to avoid padding a header-only first fragment.
- Progressive/layered encoding: earlier fragments contain earlier JPEG data for partial decoding.

Minimum fragment length:
- The smallest fragment must fit: per-fragment header (18 bytes) + per-block header (8 bytes) + 1 block (block_size bytes).
- With default block_size=223: 18 + 8 + 223 = 249 bytes => 996 bases minimum.
- 1 byte maps to 4 bases, so minimum -l is 996 bases (or higher depending on block_size).

Usage examples:
  python jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py input.jpg -o out.fasta -l 1000 --no_pad_header
  python jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py input.jpg -o out.fasta -l 1000 --no_pad_header --whiten
  python jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py input.jpg -o out.fasta -l 1200 --block_size 200
  python jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py input.jpg -o out.fasta -l 1400 --whiten --block_size 250

"""
import argparse
import math
import struct
import sys
import hashlib
import json
from datetime import datetime

# File format constants (encoder "version" of this pipeline - no RS)
MAGIC = b"DNJ6"  # Different magic to distinguish from RS version
VERSION = 6

# Global header (written once, preferably in the first fragment):
# Layout string uses big-endian ('>') and fixed-size fields:
#   4s : magic tag, e.g., b"DNJ6" (4 bytes)
#    B : version (1 byte)
#    B : block_size (1 byte) -- size of each data block in bytes
#    x : reserved (1 byte) -- was RS k, now unused
#    I : jpeg_len (4 bytes) -- unencoded JPEG length in bytes
#  32s : sha256 of original JPEG (32 bytes)
#    I : num_fragments (4 bytes)
#    I : fragment_len (4 bytes) in bases (your -l)
#    B : map_period (1 byte) -- rotating base-map period (quads)
#    B : map_seed (1 byte) -- starting phase for fragment 0
#    B : whitening enabled flag (1 byte)
#    x : reserved (1 byte) -- was rs_pad_prbs, now unused
HDR_STRUCT = struct.Struct(">4s B B B I 32s I I B B B B")

# Per-block header (prepended to every data block) identifies the block:
#    I : block_index (4 bytes)
#    I : jpeg_offset (4 bytes) -- where this block's data begins in the JPEG byte stream
BLK_STRUCT = struct.Struct(">I I")

# Per-fragment header (prepended to each fragment payload):
#    I : fragment_index (4 bytes) -- 0..N-1; -1 (0xFFFFFFFF) for header-only
#    H : num_blocks_in_fragment (2 bytes)
#    I : byte_len_of_payload (4 bytes) -- number of bytes after this header in this fragment
#   8s : nonce (8 bytes) -- used to generate per-fragment keystream
FRAG_STRUCT = struct.Struct(">I H I 8s")

# Rotating 2-bit -> base maps. We rotate which map we use as we walk across the quaternary digits.
ROT_TABLE = [
    ['A', 'C', 'G', 'T'],
    ['C', 'G', 'T', 'A'],
    ['G', 'T', 'A', 'C'],
    ['T', 'A', 'C', 'G'],
]


def chunk_bytes(data, size):
    """
    Split a bytes object into chunks of at most 'size'.
    """
    return [data[i:i+size] for i in range(0, len(data), size)]


def prng_keystream(nonce8, length):
    """
    Create a deterministic keystream (length bytes) derived from a nonce.
    This is NOT cryptographic; it's only to decorrelate structure for synthesis friendliness.

    We compute SHA-256(nonce || counter) for counter=0,1,... and concatenate until length is reached.
    """
    out = bytearray()
    counter = 0
    while len(out) < length:
        block = hashlib.sha256(nonce8 + struct.pack(">I", counter)).digest()
        out.extend(block)
        counter += 1
    return bytes(out[:length])


def xor_bytes(a, b):
    """
    XOR two bytes objects of equal length.
    """
    return bytes(x ^ y for x, y in zip(a, b))


def blocks_with_offsets(jpeg_bytes, block_size):
    """
    Slice the JPEG into block_size-byte parts and tag with original offset.
    Returns a list of tuples: (block_index, jpeg_offset, block_bytes)
    """
    parts = chunk_bytes(jpeg_bytes, block_size)
    blocks = []
    for idx, payload in enumerate(parts):
        src_off = idx * block_size
        blocks.append((idx, src_off, payload))
    return blocks


def prioritized_pack(blocks, fragment_capacity_bytes, mk_nonce_for_fragment):
    """
    Sort blocks by source JPEG offset (ascending) so earlier fragments contain earlier JPEG data.
    Greedily pack (per-block-header + block_data) units into a fragment until capacity is reached.
    For each fragment, create a per-fragment header with a per-fragment nonce.

    Returns:
      list of (fragment_index, nonce8, fragment_bytes)
      where fragment_bytes = FRAG_HEADER + concatenated block units for that fragment.
    """
    blocks_sorted = sorted(blocks, key=lambda x: x[1])
    fragments = []
    current_payload = bytearray()
    current_meta = []
    frag_index = 0
    for idx, src_off, block_data in blocks_sorted:
        unit = BLK_STRUCT.pack(idx, src_off) + block_data
        if len(current_payload) + len(unit) > fragment_capacity_bytes and len(current_payload) > 0:
            nonce8 = mk_nonce_for_fragment(frag_index)
            frag_hdr = FRAG_STRUCT.pack(frag_index, len(current_meta), len(current_payload), nonce8)
            fragments.append((frag_index, nonce8, frag_hdr + current_payload))
            frag_index += 1
            current_payload = bytearray()
            current_meta = []
        current_payload += unit
        current_meta.append((idx, src_off))
    if len(current_payload) > 0:
        nonce8 = mk_nonce_for_fragment(frag_index)
        frag_hdr = FRAG_STRUCT.pack(frag_index, len(current_meta), len(current_payload), nonce8)
        fragments.append((frag_index, nonce8, frag_hdr + current_payload))
    return fragments


def bytes_to_quads(data):
    """
    Expand each byte into four 2-bit digits (most significant first).
    Each 2-bit digit is an integer 0..3.
    """
    quads = []
    for b in data:
        quads.append((b >> 6) & 0x03)
        quads.append((b >> 4) & 0x03)
        quads.append((b >> 2) & 0x03)
        quads.append(b & 0x03)
    return quads


def quads_to_dna_rotating(quads, map_period=4, map_seed=0):
    """
    Map 2-bit digits to bases using a rotating map to reduce repeats.

    Arguments:
      quads: list of 2-bit digits (0..3)
      map_period: how many quads until the rotation cycles
      map_seed: starting rotation index for the first quad in this fragment

    Returns:
      DNA string.
    """
    dna = []
    period = max(1, map_period)
    for i, q in enumerate(quads):
        rot_index = (map_seed + i) % period
        rot_index = rot_index % 4
        dna.append(ROT_TABLE[rot_index][q])
    return "".join(dna)


def bytes_to_dna_rotating(data, map_period=4, map_seed=0):
    """
    Convenience wrapper to map bytes directly to rotating DNA bases.
    """
    return quads_to_dna_rotating(bytes_to_quads(data), map_period=map_period, map_seed=map_seed)


def estimate_fragment_capacity_bytes(length_bases):
    """
    Each byte produces 4 bases. So capacity in bytes is floor(length_bases / 4).
    """
    return length_bases // 4


def bases_from_keystream(nonce8, length, map_period=4, map_seed=0):
    """
    Create non-periodic, balanced padding bases using the same rotating mapper over a keystream.
    """
    byte_len = (length + 3) // 4
    ks = prng_keystream(nonce8, byte_len)
    pads = bytes_to_dna_rotating(ks, map_period=map_period, map_seed=map_seed)
    return pads[:length]


def main():
    parser = argparse.ArgumentParser(description="Encode JPEG into DNA WITHOUT Reed-Solomon error correction. Uses whitening, rotating mapping, and keystream padding.")
    parser.add_argument("input_jpeg", help="Path to input JPEG file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA path")
    parser.add_argument("-l", "--length", type=int, required=True, help="DNA sequence length per fragment (bases)")
    parser.add_argument("--block_size", type=int, default=223, help="Size of each data block in bytes (default 223)")
    parser.add_argument("--name_prefix", type=str, default="DNJx", help="FASTA record name prefix")
    # Rotating mapper controls
    parser.add_argument("--map_period", type=int, default=4, help="Rotation period (quads) (default 4)")
    parser.add_argument("--map_seed", type=int, default=0, help="Base mapping seed for fragment 0")
    # Whitening and padding options
    parser.add_argument("--whiten", action="store_true", help="Whiten fragment body by XOR with per-fragment keystream")
    parser.add_argument("--no_pad_header", action="store_true", help="Do not pad a header-only first fragment")
    parser.add_argument("--config", type=str, default=None, help="Optional: save encoding config to JSON file")
    args = parser.parse_args()

    if args.length <= 0:
        print("Error: --length must be positive", file=sys.stderr)
        sys.exit(1)
    if args.block_size <= 0:
        print("Error: --block_size must be positive", file=sys.stderr)
        sys.exit(1)
    if args.map_seed < 0:
        print("Error: --map_seed must be >= 0", file=sys.stderr)
        sys.exit(1)
    if args.map_period <= 0:
        print("Error: --map_period must be >= 1", file=sys.stderr)
        sys.exit(1)

    # Read JPEG
    try:
        with open(args.input_jpeg, "rb") as f:
            jpeg = f.read()
    except Exception as e:
        print("Failed to read JPEG:", e, file=sys.stderr)
        sys.exit(1)

    # Root hash seeds multiple deterministic choices (nonces, PRBS)
    root_hash = hashlib.sha256(jpeg).digest()

    # Blocks with original offsets (for priority packing) - NO RS encoding
    blocks = blocks_with_offsets(jpeg, block_size=args.block_size)

    # Capacity per fragment in bytes
    frag_capacity_bytes = estimate_fragment_capacity_bytes(args.length)

    # Minimum fragment must fit: FRAG header (18) + block header (8) + 1 block (block_size bytes)
    min_bytes = FRAG_STRUCT.size + BLK_STRUCT.size + args.block_size
    if frag_capacity_bytes < min_bytes:
        need_bases = (min_bytes * 4)
        print("Error: -l too small. Need at least {} bases (for block_size={}).".format(need_bases, args.block_size), file=sys.stderr)
        sys.exit(1)

    # Deterministic per-fragment nonce (8 bytes) derived from root hash and fragment index
    def mk_nonce_for_fragment(frag_index):
        return hashlib.sha256(root_hash + struct.pack(">I", frag_index)).digest()[:8]

    # Pack prioritized fragments
    fragments = prioritized_pack(blocks, fragment_capacity_bytes=frag_capacity_bytes, mk_nonce_for_fragment=mk_nonce_for_fragment)

    # Build global header (no RS parameters)
    num_fragments = len(fragments)
    global_header = HDR_STRUCT.pack(
        MAGIC, VERSION, args.block_size, 0,  # 0 is reserved (was RS k)
        len(jpeg), root_hash, num_fragments, args.length,
        args.map_period, args.map_seed,
        1 if args.whiten else 0, 0  # 0 is reserved (was rs_pad_prbs)
    )

    # Try to put header into fragment 0. If it doesn't fit, create a dedicated header-only fragment.
    first_frag_index, first_nonce, first_payload = fragments[0]
    if len(global_header) + len(first_payload) <= frag_capacity_bytes:
        fragments[0] = (first_frag_index, first_nonce, global_header + first_payload)
        header_only_fragment = False
    else:
        header_only_fragment = True
        header_nonce = mk_nonce_for_fragment(0xFFFFFFFF)
        header_frag_hdr = FRAG_STRUCT.pack(0xFFFFFFFF, 0, len(global_header), header_nonce)
        header_fragment_payload = header_frag_hdr + global_header
        fragments = [(-1, header_nonce, header_fragment_payload)] + fragments
        num_fragments = len(fragments)

    # Convert fragments to DNA with whitening and rotating mapping.
    # Whitening applies only to the fragment body (after FRAG_STRUCT), leaving headers legible to decoders.
    frag_dnas = []
    for i, (frag_index, frag_nonce, frag_bytes) in enumerate(fragments):
        frag_hdr = frag_bytes[:FRAG_STRUCT.size]
        frag_body = frag_bytes[FRAG_STRUCT.size:]

        if args.whiten and len(frag_body) > 0:
            ks = prng_keystream(frag_nonce, len(frag_body))
            frag_body = xor_bytes(frag_body, ks)

        whitened_payload = frag_hdr + frag_body

        # Per-fragment mapping phase: base seed + fragment index (breaks cross-fragment phase alignment)
        phase_seed = (args.map_seed + i)
        dna = bytes_to_dna_rotating(whitened_payload, map_period=args.map_period, map_seed=phase_seed)

        # Padding policy:
        if header_only_fragment and frag_index == -1 and args.no_pad_header:
            # Header-only fragment is left as-is (variable length DNA)
            frag_dnas.append(dna)
        else:
            if len(dna) < args.length:
                # Non-periodic padding derived from keystream; depends on fragment nonce and phase
                pad = bases_from_keystream(frag_nonce, args.length - len(dna),
                                           map_period=args.map_period,
                                           map_seed=phase_seed)
                dna = dna + pad
            elif len(dna) > args.length:
                dna = dna[:args.length]
            frag_dnas.append(dna)

    # Write single FASTA file with all fragments/strands (each strand is a layer)
    width = int(math.log10(max(1, len(frag_dnas)))) + 1
    try:
        with open(args.output, "w", newline="\n") as f:
            for i, dna in enumerate(frag_dnas):
                f.write(">{}_{}\n".format(args.name_prefix, str(i+1).zfill(width)))
                f.write(dna + "\n")
    except Exception as e:
        print("Failed to write FASTA:", e, file=sys.stderr)
        sys.exit(1)

    # Save config file if requested
    if args.config:
        config = {
            'date': str(datetime.now()),
            'input_jpeg': args.input_jpeg,
            'output_fasta': args.output,
            'jpeg_len': len(jpeg),
            'jpeg_sha256': root_hash.hex(),
            'num_fragments': num_fragments,
            'fragment_len_bases': args.length,
            'block_size': args.block_size,
            'map_period': args.map_period,
            'map_seed': args.map_seed,
            'whiten': args.whiten,
            'no_pad_header': args.no_pad_header,
            'header_only_fragment': header_only_fragment,
            'version': VERSION
        }
        with open(args.config, "w") as f:
            json.dump(config, f, indent=2)
        print("Config saved to", args.config)

    if header_only_fragment:
        print("Note: First fragment is header-only {}."
              .format("and unpadded" if args.no_pad_header else "(padded to -l)"))
    print("Wrote {} fragments to {}".format(len(frag_dnas), args.output))


if __name__ == "__main__":
    main()

