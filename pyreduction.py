import os
import glob
import numpy as np
import logging
import astropy.io.fits as pyfits
from pathlib import Path
from astropy.io import fits

from astropy.nddata import CCDData
from matplotlib import pyplot as plt


class Reduction:

    def __init__(self, data_dir):
        self.logger = logging.getLogger(__name__)  # Create the logger for the file
        self.data_dir = Path(data_dir)  # αυτό πρέπει να αλλάξει και να λαμβάνεται απο κάποιο yaml file απο το ui;;
        # CHECKS IF OUTPUT DIRECTORY FOR THE REDUCED FRAMES EXISTS AND CREATES IT
        self.outdir = os.path.join(self.data_dir, "reduced-frames")
        os.makedirs(self.outdir, exist_ok=True)

        # INITIALIZING FRAMES TO BE REDUCED
        # self.files = os.listdir(data_dir)
        self.bias_frames = []
        self.master_bias = []
        self.dark_frames = []
        self.flat_frames = []
        self.bias_files = []
        self.tempbias_frames = []
        self.headers = []
        self.first_frame = []
        # self.n = 0

    def get_bias(self):
        """
        Loads the bias frames to a 3D array
        """
        self.bias_files = glob.glob(os.path.join(self.data_dir, "bias*.fit"))  # Load the bias files
        self.n = len(self.bias_files)  # Get the number of bias files loaded
        self.first_bframe = fits.getdata(self.bias_files[0])  # Get the data from the first fits to initiate the 3D array
        self.imsize_y, self.imsize_x = self.first_bframe.shape  # Specify the size of the image
        self.bias_stack = np.zeros((self.imsize_y, self.imsize_x, self.n), dtype=np.float32)  # Create a 3D array of 0s
        for ff in range(0, self.n):
            self.tempbias_frames = fits.getdata(self.bias_files[ff])  # Loop through the bias list
            self.bias_stack[:, :, ff] = self.tempbias_frames  # Insert the data to the array
            # self.bias_frames.append(np.ones_like(self.tempbias_frames[0].data) * self.tempbias_frames[0].data) # Multiply an array filled with 1s to get the dimensions by the values of the fits)

    def get_master_bias(self):
        """
        Combines Bias frames with MEDIAN and writes the Master Bias
        """
        self.master_bias = np.median(self.bias_stack, axis=2)  # Calculate the median along the 2nd axis (the 3rd dimension)
        print("Median Bias: ", np.median(self.master_bias))  # print the median of the master bias
        # self.master_bias_vis = self.master_bias.reshape(890, 1024)   # (NAXIS2, NAXIS1)
        try:
            pyfits.writeto(Path(os.path.join(self.outdir, "master_bias.fit")), self.master_bias,overwrite=True)  # Write the file to the outdir path
        except:
            self.logger.exception("Something happened while writing masterbias.")

    def get_dark(self):
        """
        Loads the master bias subtracted dark frames divided by exposure time
        """
        self.dark_files = glob.glob(os.path.join(self.data_dir, "dark*.fit"))  # Load the bias files
        self.dhd = fits.getheader(self.dark_files[1])
        self.exptime = self.dhd['EXPTIME']
        self.n2 = len(self.dark_files)  # Get the number of dark files loaded
        self.first_dframe = fits.getdata(self.dark_files[0])  # Get the data from the first fits to initiate the 3D array
        self.imsize_y, self.imsize_x = self.first_dframe.shape  # Specify the size of the image
        self.dark_stack = np.zeros((self.imsize_y, self.imsize_x, self.n2), dtype=np.float32)  # Create a 3D array of 0s
        for ff in range(0, self.n2):
            self.tempdark_frames = fits.getdata(self.dark_files[ff])  # Loop through the dark list
            self.dark_stack[:, :, ff] = (self.tempdark_frames - self.master_bias) / self.exptime  # Insert the data to the array

        # self.dark_files = glob.glob(os.path.join(self.data_dir, "dark*.fit"))
        # self.dhd = fits.getheader(self.dark_files[1])
        # self.exptime = self.dhd['EXPTIME']
        # # print(self.exptime)  # Print the exposure time of the dark file for validation
        # for ff in self.dark_files:
        #     self.tempdark_frames = fits.getdata(ff)
        #     self.dark_frames.append((self.tempdark_frames - self.master_bias) / self.exptime)
        # self.dark_frames = np.array(self.dark_frames)

    def get_master_dark(self):
        """
        Combines the master bias subtracted dark frames with MEDIAN and writes the Master Dark
        """
        self.master_dark = np.median(self.dark_stack, axis=2)
        # print(self.master_dark.shape)
        for i in range(len(self.master_dark)):  # Makes sure no negative value pixels exist
            for j in range(len(self.master_dark[i])):
                if self.master_dark[i][j] < 0:
                    self.master_dark[i][j] = 0

        print("Median Dark: ", np.median(self.master_dark))  # Prints the median of the master dark
        try:  # ΕΔΩ ΝΑ ΚΟΙΤΑΞΩ ΠΙΘΑΝΟ ERROR ΠΟΥ ΘΑ ΠΡΕΠΕΙ ΝΑ ΠΕΤΑΕΙ
            pyfits.writeto(Path(os.path.join(self.outdir, "master_dark.fit")), self.master_dark, overwrite=True)
        except:
            self.logger.exception("Something happened while writing masterdark.")

    def get_flat(self):
        """
        Loads the callibrated flat frames
        """
        self.flat_files = glob.glob(os.path.join(self.data_dir, "flat*.fit"))  # Load the flat files
        print(self.flat_files)
        self.dhd2 = fits.getheader(self.flat_files[1])
        self.exptime2 = self.dhd2['EXPTIME']
        self.n3 = len(self.flat_files)  # Get the number of flat files loaded
        self.first_fframe = fits.getdata(self.flat_files[0])  # Get the data from the first fits to initiate the 3D array
        self.imsize_y, self.imsize_x = self.first_fframe.shape  # Specify the size of the image
        self.flat_stack = np.zeros((self.imsize_y, self.imsize_x, self.n3), dtype=np.float32)  # Create a 3D array of 0s
        for ff in range(0, self.n3):
            self.tempflat_frames = fits.getdata(self.flat_files[ff])  # Loop through the dark list
            self.flat_stack[:, :, ff] = self.tempflat_frames - self.master_bias - (self.exptime2 * self.master_dark)  # Insert the data to the array
        # print(self.tempflat_frames)
        # print(self.flat_stack[0, 0, 4])

        # self.flat_files = glob.glob(os.path.join(self.data_dir, "flat*.fit"))
        # self.dhd2 = fits.getheader(self.flat_files[1])
        # self.exptime2 = self.dhd2['EXPTIME']
        # for ff in self.flat_files:
        #     self.tempflat_frames = fits.getdata(ff)
        #     self.flat_frames.append(self.tempflat_frames - self.master_bias - (self.exptime2 * self.master_dark))
        # self.flat_frames = np.array(self.flat_frames)

    def get_master_flat(self):
        self.master_flat = np.median(self.flat_stack, axis=2)
        self.meanvalue = np.median(self.master_flat)
        self.master_flat_norm = self.master_flat / self.meanvalue
        print('shpae:', self.master_flat_norm.shape)
        # self.master_flat_vis = self.master_flat_norm.reshape(1024, 890)
        print(self.master_flat)
        print(self.meanvalue)
        print("norm")
        print(self.master_flat_norm)
        print("vis")
        print(self.master_flat_vis)
        try:
            pyfits.writeto(Path(os.path.join(self.outdir, "master_flat.fit")), self.master_flat, overwrite=True)
            pyfits.writeto(Path(os.path.join(self.outdir, "master_flat_norm.fit")), self.master_flat_norm, overwrite=True)
        except:
            self.logger.exception("Something happened while writing masterflat.")


data_dir = (r'/home/ggrivas/Desktop/hatp32')
test = Reduction(data_dir)
test2 = test.get_bias()
test3 = test.get_master_bias()
test4 = test.get_dark()
test5 = test.get_master_dark()
test6 = test.get_flat()
test7 = test.get_master_flat()
