#!/Users/roy/anaconda3/bin/python

# This script is used to apply NL-Means denoising on a image series
# Input: <TIF-file image series>
# Output: <TIF-file image series>

import sys
import numpy as np
from skimage import img_as_float, io
from skimage.restoration import denoise_nl_means, estimate_sigma
import multiprocessing as mp

# Read input image file
def readImg(filename):
    img = img_as_float(io.imread(filename))
    return img

# Denoise image
def denoiseImg(img, frame, patch_size=5, patch_distance=6):

    # Estimate noise standard deviation
    sigma_est = np.mean(estimate_sigma(img, channel_axis=None))
    # Apply NL-Means denoising
    denoise = denoise_nl_means(img, h=1.15 * sigma_est, fast_mode=False, patch_size=patch_size, patch_distance=patch_distance)

    # report progress
    print(f'Frame: {frame+1:>5} denoised!')

    return denoise

# Main function
def main():

    px_bitsize = 16                         # pixel bit size
    patch_size = 9                          # patch size for NL-Means
    patch_distance = 6                      # patch distance for NL-Means

    filename = sys.argv[1]                  # image file name
    img = readImg(filename)                 # read image data
    
    # store processes for parallelization
    processes = [(img[t], t, patch_size, patch_distance) for t in range(len(img))]

    # start multiprocessing
    pool = mp.Pool(mp.cpu_count())
    results = pool.starmap(denoiseImg, processes)
    pool.close

    # denoised image
    cleaned_img = np.array(results)
    
    # save output file
    outname = f"{'.'.join(filename.split('.')[:-1])}_cleaned.tif"
    io.imsave(outname, np.uint16(cleaned_img * (2**px_bitsize - 1)), check_contrast=False)

# Run main
if __name__ == '__main__':
    main()

# Ankit Roy
# 22nd November, 2022
# 23rd November, 2022
#       -->     Removing CLAHE function.
#       -->     Adding progress bar.
# 24th November, 2022
#       -->     Now runs parallel. Each frame is processed by a different CPU.
#       -->     Removed progress bar, instead reported frames denoised.
#       -->     Now does not check contrast while saving TIF files.