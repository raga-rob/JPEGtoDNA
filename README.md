### JPEGtoDNA (and back) ###

JPEGtoDNA uses a XORing and Rotatory Encoding approach to mitigate DNA synthesis contraints dure to repeats encoding repetitive DNA encoded bits of 0s and 1s.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
>bmp_to_jpeg.py: A preliminary code that converts BMP images to JPEG in case the source is BMP.

>jpeg_to_dna_rotatoryEncoding.py: The encoder is able to pick the original JPEG file and encode DNA sequences into an output .fasta file with DNA sequences that can generate an identical JPEG image to the input.

>dna_to_jpeg.py: convertes the .fasta DNA sequences back into the original JPEG image (accomodating substitutive mutations but not indels).

>dna_to_jpeg_rotatoryEncoding_progressive_layered.py: The rotatory, progressive layered encoding wit Reed-Solomon (RS) error correction is able to pick the original JPEG file and encode DNA sequences into an output .fasta file with DNA sequences that correspond to JPEG image layers of ever greater resolution to maximise the speed of meaningful data retrieval during DNA sequencing.

>jpeg_to_dna_rotatoryEncoding_progressive_layered.py: The decoder is able to pick the encoded DNA sequences in the .fasta file and produce a cumulate of layers of ever greater resolution that can generate several layers of JPEG images of ever growing resolution that builds up into the final identical JPEG image to the input when all DNA strands are decoded. This approach maximises the speed at which meaningful information retreival of the JPEG image during DNA sequencing.

> WITHOUT_RS/ : in this folder the same encoder and decoder for rotatory, progressive layered encoding is available without Reed-Solomon error correction, which saves a small amount of bytes needed to encode the image

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code is able to afford subsitutions but is not sutiable for large deletions.

Overall this code is able to encode a BMP image into a JPEG file that has many less bits and includes the progressive encoding (bmp_to_jpeg.py). Then this JPEG file can be encoded in DNA with some DNA synthesis mitigation functions such as rotatory encoding and whitening to avoid repeats and homopolymeric DNA sequences that would be difficult to make (jpeg_to_dna_rotatoryEncoding.py / dna_to_jpeg_rotatoryEncoding_progressive_layered.py / dna_to_jpeg_rotatoryEncoding_progressive_layered_no_rs.py). 

Finally when using the code (dna_to_jpeg.py / jpeg_to_dna_rotatoryEncoding_progressive_layered.py / jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py) it returns the correct JPEG image either in full (dna_to_jpeg.py) or generating images with progressive layers of greater resolution (jpeg_to_dna_rotatoryEncoding_progressive_layered.py / jpeg_to_dna_rotatoryEncoding_progressive_layered_no_rs.py).
 
Improvements yet to make with this code is that although it sustains some mutations, it cannot deal with deletions or insertions. This is okay in the context of DNA data storage because the DNA sequences can be corrected by majority vote during sequencing for multiple instances of the same sequence will be sequenced. Therefore the encoder jpeg_to_dna_rotatoryEncoding.py is already complete with the possibility to include error correction at the DNA strand level which would mitigate for potential deletions/insertions which are nevertheless unlikely due to the fact that multiple DNA strands will be available to create a consensus sequence during sequencing.
