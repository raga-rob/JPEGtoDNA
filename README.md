### JPEGtoDNA (and back) ###

JPEGtoDNA uses a XORing and Rotatory Encoding approach to mitigate DNA synthesis contraints dure to repeats encoding repetitive DNA encoded bits of 0s and 1s.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
>bmp_to_jpeg.py: A preliminary code that converts BMP images to JPEG in case the source is BMP.

>jpeg_to_dna_rotatoryEncoding.py: The decoder is able to pick the encoded DNA sequences and produce an output file as JPEG identical to the input.

>dna_to_jpeg.py: convertes the DNA sequences back into the original JPEG image (accomodating substitutive mutations but not indels).

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code is able to afford subsitutions but is not sutiable for large deletions.

Overall this code is able to encode a BMP image into a JPEG file that has many less bits and includes the progressive encoding (bmp_to_jpeg.py). Then this JPEG file can be encoded in DNA with some DNA synthesis mitigation functions such as rotatory encoding and whitening to avoid repeats and homopolymeric DNA sequences that would be difficult to make (jpeg_to_dna_rotatoryEncoding.py). 

Finally when using the code (dna_to_jpeg.py) it returns the correct JPEG image.
 
 
Improvements yet to make with this code is that although it sustains some mutations, it cannot deal with deletions or insertions. This is okay in the context of DNA data storage because the DNA sequences can be corrected by majority vote during sequencing for multiple instances of the same sequence will be sequenced. Therefore the encoder jpeg_to_dna_rotatoryEncoding.py is already complete with the possibility to include error correction at the DNA strand level which would mitigate for potential deletions/insertions which are nevertheless unlikely due to the fact that multiple DNA strands will be available to create a consensus sequence during sequencing.
