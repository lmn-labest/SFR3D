import sys

def get_newline(lines):
    global k
    k = k + 1
    return lines[k-1]

def readFile(fileIn):

    global k
    k = 0
    
    try:
        with open(fileIn,"r") as f:
            data = f.read()
    except IOError as FileNotFoundError:
        print("File {0} not found!".format(fileIn))
        sys.exit(2)

    lines = data.split('\n')

    pO2 = float(get_newline(lines).split()[1])
    pN2 = float(get_newline(lines).split()[1])
    arC = pO2,pN2
 

    name, c, h, o, n = get_newline(lines).split()
    chon = int(c), int(h), int(o), int(n)


 
    return chon,arC
