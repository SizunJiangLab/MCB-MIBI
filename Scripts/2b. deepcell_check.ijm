// yeoyaoyu, 2022

// This macro automates quality assessment of DeepCell's output

// select folder containing DeepCell output
#@ File (label = "Input directory", style = "directory") in_folder

in_folder = in_folder + "/"

file_list = newArray('nucleus.tif', 'surface.tif', 'cellOutline.tif')

for (i = 0; i < file_list.length; i++) {
	//open .tif image
	file_name = file_list[i];
	files = in_folder + file_name;
	open(files);
	run("16-bit");
}

run("Merge Channels...", "c2=surface.tif c3=nucleus.tif c4=cellOutline.tif create");
