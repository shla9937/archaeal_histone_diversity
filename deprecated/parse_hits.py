import sys


def main():
    hits_file = sys.argv[1]
    hits = import_file(hits_file)
    print('Number of hits with E <= 10:')
    print(len(hits))
    write_file(hits, hits_file)
    return True

def import_file(hits_file):
    hits = []
    ready = False
    collect = False
    f = open(hits_file, 'rt')
    for l in f:
        if 'inclusion threshold' in l:
            continue
        elif collect == True and len(l.rstrip()) == 0:
            break 
        elif collect == True:
            current_line = l.rstrip().split()
            if '|' not in current_line[8]:
                hits.append(current_line[8])
        elif 'E-value' in l:
            ready = True
        elif ready == True:
            collect = True
    return hits

def write_file(hits, hits_file):
    new_file = hits_file.rstrip('.out')+'.txt'
    f = open(new_file, 'a')
    for hit in hits:
        f.write(hit+'\n')
    return True


if __name__ == '__main__':
    main()

