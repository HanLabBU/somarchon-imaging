# Hua-an Tseng, huaantseng@gmail.com
# inspirited by https://github.com/kr-hansen/ptmc
#
# Publication:
# Population imaging of neural activity in awake behaving mice
# Kiryl D. Piatkevich, Seth Bensussen, Hua-an Tseng, Sanaya N. Shroff, Violeta Gisselle Lopez-Huerta, Demian Park, Erica E. Jung Or A. Shemesh, Christoph Straub, Howard J. Gritton, Michael F. Romano, Emma Costa, Bernardo L. Sabatini, Zhanyan Fu, Edward S. Boyden, Xue Han

import collections
import h5py
import numpy as np
import os
from scipy import ndimage
from scipy import signal
from scipy.io import loadmat
from skimage.restoration import denoise_tv_chambolle
import sys
import tifffile
import tkinter as tk
from tkinter import filedialog
from tqdm import tqdm

def motion_correction(filename_list=None, save_filename=None, save_foldername='motion_corrected', as_group = True, save_format='hdf5', **kwargs):
    if filename_list is None:
        root = tk.Tk()
        root.withdraw()
        filename_list = tk.filedialog.askopenfilenames()
        filename_list = list(filename_list)

    if len(filename_list)==1:
        print("\033[1;32;40mProcessing %s\033[0m" % filename_list[0])
    else:
        current_foldername = os.path.dirname(filename_list[0])
        print("\033[1;32;40mProcessing files in %s\033[0m" % current_foldername)

    for filename_idx, filename in enumerate(tqdm(filename_list, desc='Files')):
        [_, file_extension] = os.path.splitext(filename)
        if file_extension=='.tif':
            images = load_tif_tifffile(filename)
        elif file_extension=='.hdf5':
            file_data = h5py.File(filename, 'r')
            images = file_data['image_data']
        elif file_extension=='.mat':
            file_data = loadmat(filename)
            images = file_data['image_data']
        else:
            print("Incorrect file format.")
            break

        if as_group:
            if filename_idx==0:
                print("Generating reference with: %s" % filename)
                [image_shift, _] = calculate_shift(images, **kwargs)
                first_shifted_images = apply_shift(images, image_shift)
                ref_image = np.mean(first_shifted_images, axis=0)
            [image_shift, _] = calculate_shift(images, ref_image=ref_image, **kwargs)
        else:
            print("Generating reference with: %s" % filename)
            [image_shift, _] = calculate_shift(images, **kwargs)
            first_shifted_images = apply_shift(images, image_shift)
            ref_image = np.mean(first_shifted_images, axis=0)
            [image_shift, _] = calculate_shift(images, ref_image=ref_image, **kwargs)

        shifted_images = apply_shift(images, image_shift)

        current_foldername = os.path.dirname(filename)

        if save_filename is None:
            current_filename = os.path.basename(filename)
            current_save_filename = current_foldername+'/'+save_foldername+'/m_'+current_filename
        else:
            current_save_filename = current_foldername+'/'+save_foldername+'/'+save_filename

        [current_save_filename, _] = os.path.splitext(current_save_filename)
        current_save_filename = current_save_filename+'.'+save_format

        if os.path.isdir(current_foldername + '/' + save_foldername) is False:
            os.makedirs(current_foldername + '/' + save_foldername)

        if save_format=='npy':
            np.save(current_save_filename, shifted_images, allow_pickle=False)
        elif save_format == 'hdf5':
            opened_save_file = h5py.File(current_save_filename, "w")
            opened_save_file.create_dataset('image_data', data=shifted_images, chunks=True)
            opened_save_file.close()
        else:
            save_tif_tifffile(shifted_images, current_save_filename)

def apply_highpass_flt(images, sigma=50):
    if images.ndim == 2:
        images = np.expand_dims(images, axis=0)
        squeeze_image = True
    else:
        squeeze_image = False
    images = images.astype('float')
    # print("Applying sharpen....")
    for frame_idx, image in enumerate(images):
        lowpass_image = ndimage.gaussian_filter(image, sigma)
        images[frame_idx, :, :] = image - lowpass_image
    if squeeze_image:
        images = np.squeeze(images)
    return images

def apply_sharpen(images, sigma=[2, 1], alpha=100):
    # http://www.scipy-lectures.org/advanced/image_processing/
    if images.ndim == 2:
        images = np.expand_dims(images, axis=0)
        squeeze_image = True
    else:
        squeeze_image = False
    images = images.astype('float')
    # print("Applying sharpen....")
    for frame_idx, image in enumerate(images):
        lowpass_image = ndimage.gaussian_filter(image, sigma[0])
        filter_lowpass_image = ndimage.gaussian_filter(lowpass_image, sigma[1])
        images[frame_idx, :, :] = lowpass_image + alpha * (lowpass_image - filter_lowpass_image)
    if squeeze_image:
        images = np.squeeze(images)
    return images


def apply_shift(images, image_shift):
    if len(images) != len(image_shift):
        raise IndexError("images and image_shift require the same length in first dimension.")
    else:
        # print("Applying shift....")
        try:
            shifted_images = np.zeros(images.shape, dtype=images.dtype)  # images as numpy array
        except:
            shifted_images = np.zeros((len(images), images[0].shape[0], images[0].shape[1]), dtype='uint16')  # images as nd2

        for frame_idx, image in enumerate(tqdm(images, desc="Applying shift", leave=False)):
            if np.array_equal(image_shift[frame_idx], [0, 0]):
                shifted_images[frame_idx, :, :] = image
            else:
                shifted_images[frame_idx, :, :] = ndimage.shift(image, image_shift[frame_idx])

        return shifted_images


def apply_std_eq(images):
    if images.ndim == 2:
        images = np.expand_dims(images, axis=0)
        squeeze_image = True
    else:
        squeeze_image = False
    image_shape = images.shape[1:]
    # print("Applying histogram equalization....")
    for frame_idx, image in enumerate(images):
        images[frame_idx, :, :] = (image - image.mean()) / image.std()
    if squeeze_image:
        images = np.squeeze(images)
    return images

def calculate_shift(images, ref_image=None, process_functions='voltage'):
    # images is numpy array frame x height x width
    # global_subpixel: don't use it
    # local_subpixel: most of time, it's not better

    if process_functions=='voltage':
        process_functions=['remove_edges', 'apply_highpass_flt', 'apply_sharpen', 'remove_intensity', 'apply_std_eq']
    else:
        pass

    process_function_list = {
        'apply_highpass_flt': apply_highpass_flt, 
        'apply_sharpen': apply_sharpen,
        'apply_std_eq': apply_std_eq,
        'remove_edges': remove_edges,
        'remove_intensity': remove_intensity,
    }

    images = np.array(images)

    if images.ndim == 2:
        images = np.expand_dims(images, axis=0)

    if ref_image is None:
        ref_image = np.mean(images, axis=0)
    elif ref_image=='first':
        ref_image = images[0]

    print(" ")
    print("Image size: %sx%s, frame numeber: %s" % (images.shape[1], images.shape[2], images.shape[0]))
    print("Process functions (%s): %s" % (len(process_functions), ', '.join(map(str, process_functions))))

    image_shift = np.zeros((len(images), 2))

    if process_functions is not None:
        for process_function_name in process_functions:
            current_process_function = process_function_list[process_function_name]
            ref_image = current_process_function(ref_image)

    if global_subpixel < 1:
        ref_image = ndimage.zoom(ref_image, (1 / subpixel, 1 / subpixel), order=1)

    ref_image_center = np.array(ref_image.shape) // 2
    ref_image_fft = np.fft.fft2(ref_image)

    # print("Calculating shift....")
    for frame_idx, image in enumerate(tqdm(images, desc="Calculating shift", leave=False)):

        if process_functions is not None:
            for process_function_name in process_functions:
                current_process_function = process_function_list[process_function_name]
                image = current_process_function(image)

        image_fft = np.fft.fft2(image)
        cross_correlation = abs(np.fft.ifft2(image_fft * ref_image_fft.conjugate()))
        cross_correlation_peak = np.array(
            np.unravel_index(np.argmax(np.fft.fftshift(cross_correlation)), ref_image.shape))
        image_shift[frame_idx, :] = -1 * (cross_correlation_peak - ref_image_center)

    return image_shift, ref_image


def load_tif_tifffile(filename_list):
    # load all frames
    # fast for frame number >1000
    # check load_tif_pil(filename_list, start_frame, frame_number)
    if isinstance(filename_list, str):
        filename_list = [filename_list]
    output_images = collections.deque([])
    # print("Loading files....")
    for filename in filename_list:
        with tifffile.TiffFile(filename) as image_data:
            temp_images = image_data.asarray()
            output_images.append(temp_images)
    output_images = np.dstack(output_images)
    return output_images  # frame x height x width

def remove_intensity(images, low_std_threshold=-1, high_std_threshold=1):
    # remove any pixel with intensity below mean+(std_threshold*std)
    if images.ndim == 2:
        images = np.expand_dims(images, axis=0)
        squeeze_image = True
    else:
        squeeze_image = False
    images = images.astype('float')
    for frame_idx, image in enumerate(images):
        std_image = np.std(image)
        image[image<(np.mean(image)+low_std_threshold*std_image)] = 0
        image[image>(np.mean(image)+high_std_threshold*std_image)] = np.mean(image)+high_std_threshold*std_image
        images[frame_idx, :, :] = image
    if squeeze_image:
        images = np.squeeze(images)
    return images


def remove_edges(images, edges=0.1):
    # remove any pixel at the edges: UDLR
    if images.ndim == 2:
        images = np.expand_dims(images, axis=0)
        squeeze_image = True
    else:
        squeeze_image = False
    images = images.astype('float')

    edges = np.array(edges)
    if edges.size == 1:
        edges = np.repeat(edges, 4)

    if any(edges < 1):
        u_pixel = int(np.floor(images.shape[1] * edges[0]))
        d_pixel = int(np.floor(images.shape[1] * edges[1]))
        l_pixel = int(np.floor(images.shape[2] * edges[2]))
        r_pixel = int(np.floor(images.shape[2] * edges[3]))
    else:
        u_pixel = int(np.floor(edges[0]))
        d_pixel = int(np.floor(edges[1]))
        l_pixel = int(np.floor(edges[2]))
        r_pixel = int(np.floor(edges[3]))
    # print("Cropping edges: U%s D%s L%s R%s" % (u_pixel,d_pixel,l_pixel,r_pixel))

    images = images[:, u_pixel + 1:-1 * d_pixel - 1, l_pixel + 1:-1 * r_pixel - 1]

    if squeeze_image:
        images = np.squeeze(images)
    return images

def save_tif_pil(images, save_filename):
    if images.ndim == 2:
        images = PIL.Image.fromarray(images)
        images.save(save_filename, save_all=True)
    else:
        saved_image = []
        for image in images:
            saved_image.append(PIL.Image.fromarray(image))

        saved_image[0].save(save_filename, save_all=True, append_images=saved_image[1:])

def save_tif_tifffile(images, save_filename):
    tifffile.imsave(save_filename, images,photometric='minisblack')







