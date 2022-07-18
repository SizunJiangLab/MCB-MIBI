// yeoyaoyu, 2022

// This macro automates bit size conversion and pixel normalization that are required for cell segmentation.

// select input/output folder for nucleus and surface .tif
#@ File (label = "Input directory", style = "directory") in_folder

in_folder = in_folder + "/"

file_list = getFileList(in_folder);

for (i = 0; i < file_list.length; i++) {
	//open .tif image
	file_name = file_list[i];
	files = in_folder + file_name;
	open(files);	
	//set 16-bit
	run("16-bit");
	//normalize (0.3% saturation)
	run("Enhance Contrast...", "saturated=0.3 normalize");
	//save as .tif
	if (file_name == 'DAPI.tif') {
		run("Blue");
		saveAs("Tiff", in_folder + "nucleus.tif");
	} else {
		run("Green");
		saveAs("Tiff", in_folder + "surface.tif");
	}
	}
