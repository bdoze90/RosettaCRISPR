"""File to test simple functions"""


from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

def replace(file_path):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_line = line[:10] + new_chain + line[11:]
                new_file.write(new_line)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)






replace("/Users/brianmendoza/Desktop/testfile.txt")