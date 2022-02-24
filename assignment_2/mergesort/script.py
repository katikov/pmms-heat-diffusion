def calculate(path):
    count = 0
    sum = 0
    key = 0
    dic = {}
    fname = path
    fh = open(fname+'.txt')
    for line in fh:
        if line[0].isdigit() :
            key = line.strip('\n')
            count = 0
            sum = 0
        elif line[0].isalpha() :
            value = line.split(':')[-1]
            value = value.strip(' seconds \n')
            value = float(value)
            sum += value
            count += 1
            if(count==3):
                dic[key] = sum/3
    print('-----'+fname+'-----')
    for x in dic:
        print(x,':',dic[x])
    print('')

if __name__ == '__main__':
    file_path = ['result_threshold_dynamic',
                 'result_threshold_guided',
                 'result_threshold_static',
                 'result_threshold',
                 'result4096_chunk_dynamic',
                 'result4096_chunk_guided',
                 'result4096_chunk_static',
                 'result_inner',
                 'result_outer',
                 'result_merge_par',
                 'result_merge_seq',]
    for path in file_path:
        calculate(path)