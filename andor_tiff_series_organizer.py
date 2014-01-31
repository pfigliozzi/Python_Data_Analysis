import os
import string
import Tkinter,tkFileDialog

'''This script will clean up a folder that has both the andor raw '.sif' format
and the exported '.tif' files. This is will created a folder based on the names
of the raw data and organize all the cooresponding tifs into the appropriate 
folders while renaming the files such that their names are only numbers.'''

root = Tkinter.Tk()
root.withdraw()
directory = tkFileDialog.askdirectory()
for files in os.listdir(directory):
	full_file_path=os.path.join(directory, files)
	if files.endswith('.sif') and os.stat(full_file_path).st_size > 1000000:
		len_filename=len(os.path.splitext(files)[0])
		new_directory=os.path.join(directory, os.path.splitext(files)[0])
		try:
			os.mkdir(new_directory)
		finally:
			filename=os.path.splitext(files)[0]
			for file in glob.glob(os.path.join(directory,filename+'*.tif')):
				filename=os.path.split(file)[1]
				new_name='1'+filename[len_filename+1:]
				os.rename(file, os.path.join(new_directory,new_name))
				
