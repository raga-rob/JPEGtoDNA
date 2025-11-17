"""
bmp_to_jpeg.py

Converts a BMP image to a progressive JPEG with configurable quality.
Progressive JPEG places image data in multiple scans that allow early,
coarse previews that refine as more bytes are decoded.

Usage:
  python bmp_to_jpeg.py input.bmp -o output.jpg
  python bmp_to_jpeg.py input.bmp -o output.jpg --quality 90
"""

import argparse
import io
import sys
from PIL import Image


def bmp_to_progressive_jpeg_bytes(bmp_path, quality=85):
    """
    Open a BMP image from disk, convert to RGB, and encode as progressive JPEG.

    Arguments:
      bmp_path: path to the input BMP file
      quality: JPEG quality (1..95 typical). Higher is better quality and larger size.

    Returns:
      Bytes object containing the encoded JPEG file.
    """
    with Image.open(bmp_path) as img:
        # Ensure a standard three-channel RGB; BMP may carry palette/other modes.
        img = img.convert("RGB")
        out = io.BytesIO()
        # progressive=True writes a multi-scan JPEG, useful for "priority" retrieval.
        # optimize=True lets Pillow optimize Huffman tables for a slightly smaller file.
        img.save(out, format="JPEG", progressive=True, quality=quality, optimize=True)
        return out.getvalue()


def main():
    parser = argparse.ArgumentParser(description="Convert BMP to progressive JPEG.")
    parser.add_argument("input_bmp", help="Path to input BMP file")
    parser.add_argument("-o", "--output", required=True, help="Output JPEG file path")
    parser.add_argument("--quality", type=int, default=85, help="JPEG quality (default 85)")
    args = parser.parse_args()

    try:
        jpeg_bytes = bmp_to_progressive_jpeg_bytes(args.input_bmp, quality=args.quality)
    except Exception as e:
        print("Failed to convert BMP to JPEG:", e, file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.output, "wb") as f:
            f.write(jpeg_bytes)
    except Exception as e:
        print("Failed to write output JPEG:", e, file=sys.stderr)
        sys.exit(1)

    print("Wrote progressive JPEG to", args.output)


if __name__ == "__main__":
    main()
