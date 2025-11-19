"""
dna_to_jpeg_rotatoryEncoding_progressive_layered.py

Decodes DNA FASTA back into JPEG with progressive/layered decoding support.
- Supports partial decoding: generates multiple partial images as more fragments are decoded
- Each layer adds more fragments to progressively improve image quality
- Uses Reed-Solomon error correction for robustness
- Handles rotating base mapping and whitening

Usage examples:
  python dna_to_jpeg_rotatoryEncoding_progressive_layered.py input.fasta -o output.jpg
  python dna_to_jpeg_rotatoryEncoding_progressive_layered.py input.fasta -o output.jpg --partially_decoded_layers 5
  python dna_to_jpeg_rotatoryEncoding_progressive_layered.py input.fasta -o output.jpg --partially_decoded_layers 10 --partial_output_dir partial_images

Install dependencies:
  pip install reedsolo
"""

import argparse
import struct
import sys
import math
import hashlib
import os

try:
    import reedsolo
except ImportError:
    print("Missing dependency: reedsolo. Install with: pip install reedsolo", file=sys.stderr)
    raise

# Constants
MAGIC = b"DNJ5"
VERSION_SUPPORTED = {5}
HDR_STRUCT = struct.Struct(">4s B B B I 32s I I B B B B")
BLK_STRUCT = struct.Struct(">I I")
FRAG_STRUCT = struct.Struct(">I H I 8s")

ROT_TABLE = [
    ['A', 'C', 'G', 'T'],
    ['C', 'G', 'T', 'A'],
    ['G', 'T', 'A', 'C'],
    ['T', 'A', 'C', 'G'],
]


def parse_fasta(path):
    """
    Parse FASTA file and return list of (header, sequence) tuples.
    """
    records = []
    current_header = None
    current_seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, "".join(current_seq)))
                current_header = line[1:].strip()
                current_seq = []
            else:
                s = "".join(ch for ch in line if ch in "ACGT")
                current_seq.append(s)
        if current_header is not None:
            records.append((current_header, "".join(current_seq)))
    return records


def quads_from_dna_rotating(dna, map_period, map_seed):
    """
    Reverse of rotating DNA mapping: convert DNA back to quads.
    """
    rev_maps = []
    for m in ROT_TABLE:
        rev = {m[d]: d for d in range(4)}
        rev_maps.append(rev)

    period = max(1, map_period)
    quads = []
    for i, base in enumerate(dna):
        rot_index = (map_seed + i) % period
        rot_index = rot_index % 4
        q = rev_maps[rot_index][base]
        quads.append(q)
    return quads


def bytes_from_quads(quads):
    """
    Convert quads (2-bit digits) back to bytes.
    Groups every 4 quads into one byte.
    """
    out = []
    n = len(quads) // 4
    for i in range(n):
        q0 = quads[4*i + 0] & 0x03
        q1 = quads[4*i + 1] & 0x03
        q2 = quads[4*i + 2] & 0x03
        q3 = quads[4*i + 3] & 0x03
        b = (q0 << 6) | (q1 << 4) | (q2 << 2) | q3
        out.append(b)
    return bytes(out)


def prng_keystream(nonce8, length):
    """
    Create a deterministic keystream (length bytes) derived from a nonce.
    """
    out = []
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


def decode_fragment_bytes(dna, map_period, map_seed):
    """
    Decode DNA fragment to bytes using rotating mapping.
    """
    quads = quads_from_dna_rotating(dna, map_period, map_seed)
    data = bytes_from_quads(quads)
    return data


def recover_global_header(fragments_bytes, whiten_flag, map_period, map_seed_base):
    """
    Find and parse the global header, handling whitening if needed.
    """
    for list_index, frag in enumerate(fragments_bytes):
        if len(frag) < FRAG_STRUCT.size + HDR_STRUCT.size:
            continue
            
        frag_hdr = frag[:FRAG_STRUCT.size]
        frag_body = frag[FRAG_STRUCT.size:]
        
        # If whitening is enabled, we need to unwhiten the fragment body first
        if whiten_flag == 1:
            frag_index, num_blocks_in_fragment, body_len_bytes, nonce8 = FRAG_STRUCT.unpack(frag_hdr)
            if len(frag_body) >= body_len_bytes:
                # Unwhiten the fragment body
                ks = prng_keystream(nonce8, body_len_bytes)
                frag_body = xor_bytes(frag_body[:body_len_bytes], ks)
        
        # Check if body starts with MAGIC
        if len(frag_body) >= 4 and frag_body[:4] == MAGIC:
            if len(frag_body) >= HDR_STRUCT.size:
                header = HDR_STRUCT.unpack(frag_body[:HDR_STRUCT.size])
                return header, list_index
    
    raise ValueError("Global header not found in any fragment")


def decode_jpeg_from_fragments(fragments_bytes, map_period, map_seed_base, whiten_flag, rs_n, rs_k, jpeg_len):
    """
    Decode JPEG from fragment bytes.
    Returns the decoded JPEG bytes.
    """
    # Extract fragments and collect block units
    block_units = []
    rs_codeword_len = rs_n  # 255

    for list_index, frag in enumerate(fragments_bytes):
        if len(frag) < FRAG_STRUCT.size:
            continue  # Skip invalid fragments
            
        frag_index, num_blocks_in_fragment, body_len_bytes, nonce8 = FRAG_STRUCT.unpack(frag[:FRAG_STRUCT.size])
        body = frag[FRAG_STRUCT.size:]
        frag_body = body[:body_len_bytes] if len(body) >= body_len_bytes else body

        # If this fragment's body begins with MAGIC, it is the global header fragment
        if len(frag_body) >= 4 and frag_body[:4] == MAGIC:
            frag_body = frag_body[HDR_STRUCT.size:]

        # Unwhiten if required
        if whiten_flag == 1 and len(frag_body) > 0:
            ks = prng_keystream(nonce8, len(frag_body))
            frag_body = xor_bytes(frag_body, ks)

        # Parse block units
        cursor = 0
        for _ in range(num_blocks_in_fragment):
            if cursor + BLK_STRUCT.size + rs_codeword_len > len(frag_body):
                break  # Truncated fragment, skip remaining blocks
            blk_index, jpeg_offset = BLK_STRUCT.unpack(frag_body[cursor:cursor+BLK_STRUCT.size])
            cursor += BLK_STRUCT.size
            if cursor + rs_codeword_len > len(frag_body):
                break
            codeword = frag_body[cursor:cursor+rs_codeword_len]
            cursor += rs_codeword_len
            block_units.append((jpeg_offset, blk_index, codeword))

    # RS decode and assemble
    parity_symbols = rs_n - rs_k
    rs = reedsolo.RSCodec(parity_symbols)
    
    output = [0] * jpeg_len
    block_units.sort(key=lambda t: t[0])

    for jpeg_offset, blk_index, codeword in block_units:
        try:
            decoded_tuple = rs.decode(codeword)
            msg = decoded_tuple[0]
            
            if isinstance(msg, (bytes, bytearray)):
                msg_list = list(msg)
            else:
                msg_list = list(msg)
        except Exception as e:
            # RS decode failed, skip this block
            continue

        end = min(jpeg_offset + rs_k, jpeg_len)
        if end > jpeg_offset:
            output[jpeg_offset:end] = msg_list[:end - jpeg_offset]

    return bytes(output)


def decode_jpeg_from_fasta(fasta_path, out_path, verify_sha256=True):
    """
    Decode JPEG from FASTA file.
    
    Args:
        fasta_path: Path to input FASTA file
        out_path: Path to output JPEG file
        verify_sha256: Whether to verify SHA-256 hash (True for full decode, False for partial layers)
    """
    # 1) Read FASTA
    records = parse_fasta(fasta_path)
    if not records:
        raise ValueError("FASTA has no sequences")

    # 2) First pass: decode fragments with default parameters to find header
    found_header = None
    found_params = None
    
    # Try different parameter combinations
    for map_period in [4, 8, 16]:
        for map_seed_base in range(4):
            for whiten_flag in [0, 1]:
                try:
                    # Decode all fragments with these parameters
                    fragments_bytes = []
                    for list_index, (hdr, dna) in enumerate(records):
                        phase_seed = map_seed_base + list_index
                        frag_bytes = decode_fragment_bytes(dna, map_period=map_period, map_seed=phase_seed)
                        fragments_bytes.append(frag_bytes)
                    
                    # Try to find header with these parameters
                    header_tuple, header_list_idx = recover_global_header(fragments_bytes, whiten_flag, map_period, map_seed_base)
                    
                    # Verify this looks like a valid header
                    magic, version, rs_n, rs_k, jpeg_len, jpeg_sha256, num_fragments_hdr, fragment_len_bases, map_period_hdr, map_seed_base_hdr, whiten_flag_hdr, rs_pad_prbs_flag = header_tuple
                    
                    if (magic == MAGIC and 
                        version in VERSION_SUPPORTED and 
                        rs_n == 255 and
                        map_period_hdr == map_period and
                        map_seed_base_hdr == map_seed_base and
                        whiten_flag_hdr == whiten_flag):
                        found_header = header_tuple
                        found_params = (map_period, map_seed_base, whiten_flag)
                        break
                        
                except Exception:
                    continue
            
            if found_header is not None:
                break
        if found_header is not None:
            break
    
    if found_header is None:
        raise ValueError("Failed to locate valid header; ensure input matches encoder v5")

    # Parse the found header
    magic, version, rs_n, rs_k, jpeg_len, jpeg_sha256, num_fragments_hdr, fragment_len_bases, map_period, map_seed_base, whiten_flag, rs_pad_prbs_flag = found_header

    # 3) Second pass: decode all fragments with correct parameters
    fragments_bytes = []
    for list_index, (hdr, dna) in enumerate(records):
        phase_seed = map_seed_base + list_index
        frag_bytes = decode_fragment_bytes(dna, map_period=map_period, map_seed=phase_seed)
        fragments_bytes.append(frag_bytes)

    # 4) Decode with all fragments
    output_bytes = decode_jpeg_from_fragments(
        fragments_bytes, map_period, map_seed_base, whiten_flag, rs_n, rs_k, jpeg_len
    )
    
    # 5) Verify and write
    digest = hashlib.sha256(output_bytes).digest()
    if verify_sha256 and digest != jpeg_sha256:
        print("Warning: SHA-256 mismatch; reconstructed JPEG may be corrupt", file=sys.stderr)
        print(f"  Expected: {jpeg_sha256.hex()}", file=sys.stderr)
        print(f"  Got:      {digest.hex()}", file=sys.stderr)
    elif verify_sha256:
        print("SHA-256 verification passed")

    with open(out_path, "wb") as f:
        f.write(output_bytes)


def decode_jpeg_progressive_layers(fasta_path, output_dir=None):
    """
    Decode JPEG progressively - each strand in FASTA is a layer.
    Automatically decodes all layers: strand 1 = layer 1, strands 1-2 = layer 2, etc.
    
    Args:
        fasta_path: Path to input FASTA file
        output_dir: Directory for partial images (default: current directory)
    """
    # Read all records from FASTA
    records = parse_fasta(fasta_path)
    if not records:
        raise ValueError("FASTA has no sequences")
    
    total_fragments = len(records)
    
    if output_dir is None:
        output_dir = "."
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Decoding {total_fragments} progressive layers (one per strand)...")
    
    # Decode each layer progressively - each strand is a layer
    for layer in range(1, total_fragments + 1):
        # Use first N fragments for this layer (like jpeg_strand_encoding_1.py)
        layer_records = records[:layer]
        
        # Create temporary FASTA content for this layer
        layer_fasta_content = ""
        for hdr, dna in layer_records:
            layer_fasta_content += f">{hdr}\n{dna}\n"
        
        # Write to temporary file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
            tmp.write(layer_fasta_content)
            tmp_fasta = tmp.name
        
        try:
            # Output filename for this layer
            partial_filename = os.path.join(output_dir, f"partially_decoded_jpeg_{layer}.jpg")
            
            # Decode this layer (skip SHA-256 verification for partial layers)
            decode_jpeg_from_fasta(
                tmp_fasta,
                partial_filename,
                verify_sha256=False  # Partial layers won't match SHA-256
            )
            print(f"  Layer {layer}/{total_fragments}: strand {layer} -> {partial_filename}")
        except Exception as e:
            print(f"  Warning: Failed to decode layer {layer}: {e}", file=sys.stderr)
        finally:
            # Clean up temp file
            if os.path.exists(tmp_fasta):
                os.unlink(tmp_fasta)


def main():
    parser = argparse.ArgumentParser(description="Decode DNA FASTA back into the original JPEG (v5) with progressive/layered decoding. Each strand in the FASTA is automatically decoded as a separate layer.")
    parser.add_argument("input_fasta", help="Path to input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output JPEG file path (full decode from all strands)")
    parser.add_argument("--partial_output_dir", type=str, default=None,
                        help="Directory for partial layer images (default: same directory as output)")
    parser.add_argument("--no_layers", action="store_true", 
                        help="Skip progressive layer decoding, only do full decode")
    args = parser.parse_args()

    try:
        # Set up output directory for partial images
        if args.partial_output_dir is None:
            args.partial_output_dir = os.path.dirname(args.output) or "."
        
        # Decode progressive layers (each strand = one layer)
        if not args.no_layers:
            decode_jpeg_progressive_layers(
                args.input_fasta,
                args.partial_output_dir
            )
        
        # Also decode full FASTA (all strands)
        print(f"\nDecoding full FASTA (all strands): {args.input_fasta}")
        decode_jpeg_from_fasta(
            args.input_fasta,
            args.output
        )
        print(f"Full decode: {args.input_fasta} -> {args.output}")
            
    except Exception as e:
        print("Decode failed:", e, file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

