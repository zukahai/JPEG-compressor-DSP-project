# JPEG-Compression-Decompression
Implementation of the JPEG Image Compression and Decompression Algorithm
Description of the main modules:
jpegenco: It compresses an image using the JPEG Compression algorithm given the appropriate inputs.
jpegdeco: It decompresses the image that was compressed using the JPEG Compression algorithm.
Description of the sub modules:
There is no need to use the submodules, the main modules can be used directly.
myDCT: It takes an image as an input divide it into 8 * 8 blocks and compute the Discrete Cosine Transform of each block.
myIDCT: It takes a transformed image using myDCT as an input and computes the Inverse Discrete Cosine Transform of each block.
EncodeRL: Encodes a binary stream using run length encoding.
DecodeRL: Decodes a binary stream that is encoded using runlength encoding.
zigzag_pattern: converts a two dimension block into a one dimensional array that is equivalent to reading the two dimension block in a serpentine pattern which is an important part of the JPEG Operation.
jpeg_quantize: divides each 8 * 8 block by the quantization table with controls the degree of compression the rounds the result.
jpeg_dequantize: It reverses the operation of jpeg_quantize by multiplying each 8 * 8 block by the quantization table.
