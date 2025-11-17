"""
dna_to_jpeg.py - Fixed version
"""

import argparse
import struct
import sys
import math
import hashlib

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
    rev_maps = []
    for m in ROT_TABLE:
        rev = { m[d]: d for d in range(4) }
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
    out = []
    counter = 0
    while len(out) < length:
        block = hashlib.sha256(nonce8 + struct.pack(">I", counter)).digest()
        out.extend(block)
        counter += 1
    return bytes(out[:length])

def xor_bytes(a, b):
    return bytes(x ^ y for x, y in zip(a, b))

def decode_fragment_bytes(dna, map_period, map_seed):
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

def decode_jpeg_from_fasta(fasta_path, out_path):
    # 1) Read FASTA
    records = parse_fasta(fasta_path)
    if not records:
        raise ValueError("FASTA has no sequences")

    # 2) First pass: decode fragments with default parameters to find header
    # We need to try different combinations of parameters to find the header
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

    # 4) Extract fragments and collect block units
    block_units = []
    rs_codeword_len = rs_n  # 255

    for list_index, frag in enumerate(fragments_bytes):
        if len(frag) < FRAG_STRUCT.size:
            raise ValueError("Fragment {} too short for fragment header".format(list_index))
        frag_index, num_blocks_in_fragment, body_len_bytes, nonce8 = FRAG_STRUCT.unpack(frag[:FRAG_STRUCT.size])
        body = frag[FRAG_STRUCT.size:]
        frag_body = body[:body_len_bytes]

        # If this fragment's body begins with MAGIC, it is the global header fragment
        if frag_body[:4] == MAGIC:
            frag_body = frag_body[HDR_STRUCT.size:]

        # Unwhiten if required
        if whiten_flag == 1 and len(frag_body) > 0:
            ks = prng_keystream(nonce8, len(frag_body))
            frag_body = xor_bytes(frag_body, ks)

        # Parse block units
        cursor = 0
        for _ in range(num_blocks_in_fragment):
            if cursor + BLK_STRUCT.size + rs_codeword_len > len(frag_body):
                raise ValueError("Fragment {} body truncated while reading blocks".format(list_index))
            blk_index, jpeg_offset = BLK_STRUCT.unpack(frag_body[cursor:cursor+BLK_STRUCT.size])
            cursor += BLK_STRUCT.size
            codeword = frag_body[cursor:cursor+rs_codeword_len]
            cursor += rs_codeword_len
            block_units.append((jpeg_offset, blk_index, codeword))

    # 5) RS decode and assemble
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
            raise ValueError("RS decode failed for block {} at offset {}: {}".format(blk_index, jpeg_offset, e))

        end = min(jpeg_offset + rs_k, jpeg_len)
        if end > jpeg_offset:
            output[jpeg_offset:end] = msg_list[:end - jpeg_offset]

    # 6) Verify and write
    output_bytes = bytes(output)
    digest = hashlib.sha256(output_bytes).digest()
    if digest != jpeg_sha256:
        raise ValueError("SHA-256 mismatch; reconstructed JPEG may be corrupt")

    with open(out_path, "wb") as f:
        f.write(output_bytes)

def main():
    parser = argparse.ArgumentParser(description="Decode DNA FASTA back into the original JPEG (v5).")
    parser.add_argument("input_fasta", help="Path to input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output JPEG file path")
    args = parser.parse_args()

    try:
        decode_jpeg_from_fasta(args.input_fasta, args.output)
    except Exception as e:
        print("Decode failed:", e, file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print("Wrote JPEG to", args.output)

if __name__ == "__main__":
    main()
