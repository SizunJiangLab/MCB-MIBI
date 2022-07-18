// yeoyaoyu, 2022

// This macro automates the preparation of surface.tif files that are required for cell segmentation.

// select input/output folder for surface markers
#@ File (label = "Input directory", style = "directory") in_folder
#@ File (label = "Output directory", style = "directory") out_folder

in_folder = in_folder + "/"
out_folder = out_folder + "/"

// open all .tif images in the input (surface) folder
file_list = getFileList(in_folder);
for (i = 0; i < file_list.length; i++) {
	file_name = file_list[i];
	surface_file = in_folder + file_name;
	open(surface_file);	
	}

// merge and stack all .tif images
run("Images to Stack", "name=surface title=[] use");
run("Z Project...", "projection=[Sum Slices]");

// Save as .tif
saveAs("Tiff", out_folder + "merged_surface.tif")
