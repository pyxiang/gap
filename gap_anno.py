
import sys
from intervaltree import Interval, IntervalTree


def interTree(filename):   
    t = IntervalTree()
    count = 0
    with open(filename, 'r') as f :       
        for each_line in f:           
            line = each_line.strip().split('\t')
            chrname = line[1]
            posbegin = int(line[2])
            posend = int(line[3])
            gaptype = line[7]
            countname = 'gap' + str(count)
            data = [chrname, gaptype, countname]
            t[posbegin:posend] = data
            count += 1
        return(t)


class bnd(object):
    def __init__(self,bnd_string):
        self.bnd_string = bnd_string
        self.bnd_string = self.bnd_string.replace(">","")
        self.bnd_string = self.bnd_string.replace("<","")
        if '[' in self.bnd_string:
            first_right_braket = self.bnd_string.index('[')
            next_right_braket = self.bnd_string.index('[', first_right_braket+1)
            if first_right_braket != 0:
                self.stat = "s1"
            else:
                self.stat = "s4"
            self.pos = self.bnd_string[first_right_braket+1:next_right_braket]
        elif ']' in self.bnd_string:
            first_left_braket = self.bnd_string.index(']')
            next_left_braket = self.bnd_string.index(']', first_left_braket+1)
            if first_left_braket != 0:
                self.stat = "s2"
            else:
                self.stat = "s3"
            self.pos = self.bnd_string[first_left_braket+1:next_left_braket]
        self.chrom = self.pos.split(":")[0]
        self.pos_num = self.pos.split(":")[1]
      
class Getvcfinfo():
    def __init__(self, each_line):
        line = each_line.strip().split('\t')
        chrname = line[0]
        posbegin = int(float(line[1]))
        posend = line[7].split(';')[3].split('=')[1]
        chrname2 = line[7].split(';')[2].split('=')[1]
        ID = line[2]
        svtype = line[4]
        if '.' in posend:
            b = bnd(each_line)
            posend = b.pos_num
            chrname2 = b.chrom
        posend = int(float(posend))
        if 'chr' not in chrname:
            chrname = 'chr' + chrname
        if 'chr' not in chrname2:
            chrname2 = 'chr' + chrname2
        self.chrname = chrname
        self.posbegin = posbegin
        self.posend  = posend
        self.chrname2 = chrname2
        self.ID = ID            
        self.svtype = svtype

def calTemp(temp, a_chrname, a_chrname2, a_posbegin, a_posend, pos1, pos2, pos3, pos4):
    if temp:        
        temp = list(temp)
        region = ''
        datasetid = []
        gaptype = []
        for i in temp:
            if i[2][0] == a_chrname == a_chrname2:
                gaptype.append(i[2][1])
                datasetid.append(i[2][2])              
                region = '{0}:{1}-{2}'.format(i[2][0], str(i[0]), str(i[1]))
            elif i[2][0] == a_chrname:
                if i[0] <= pos1 and pos2 <= i[1]:
                    gaptype.append(i[2][1])
                    datasetid.append(i[2][2])
                    region = '{0}:{1}-{2}'.format(i[2][0], str(i[0]), str(i[1]))
            elif i[2][0] == a_chrname2:
                if i[0] <= pos3 and pos4 <= i[1]:
                    gaptype.append(i[2][1])
                    datasetid.append(i[2][2])
                    region = '{0}:{1}-{2}'.format(i[2][0], str(i[0]), str(i[1]))
                    
        return region, gaptype, datasetid
                   
    else:
        return '', '', ''
                  
def calRegion(t, a, each_line2,outfile, bp):
    region1 = ''
    region3 = ''
    datasetid1 = ''
    datasetid3 = ''
    gaptype1 = ''
    gaptype3 = ''
    if a.chrname == a.chrname2:
        bp = 1
    else:
        bp = bp
    pos1 = a.posbegin - bp
    pos2 = a.posbegin + bp
    pos3 = a.posend - bp
    pos4 = a.posend + bp
    if t.overlaps(pos1, pos2) or t.overlaps(pos3, pos4):
        temp = t.overlap(pos1, pos2)
        temp3 = t.overlap(pos3, pos4)
        region1, gaptype1, datasetid1 = calTemp(temp, a.chrname, a.chrname2, a.posbegin, a.posend, pos1,pos2,pos3, pos4)
        region3, gaptype3, datasetid3 = calTemp(temp3, a.chrname, a.chrname2, a.posbegin, a.posend, pos1, pos2,pos3, pos4)
              
    return datasetid1, datasetid3, region1, region3, gaptype1, gaptype3

def main():
 
    filename2 = sys.argv[1]
    outfile = sys.argv[2]
    bp = int(sys.argv[3])
    dataset = sys.argv[4]
    t = interTree(dataset)
    
    flag = 1
    count_vcf = 0
    count2 = 0
    datasetid2 = []
    datasetid4 = []
    gaptype2 = []
    gaptype4 = []

    with open(filename2,'r') as f2:
        for each_line2 in f2:
            if each_line2.strip()[0] == '#':
                continue
            else:
                count2 += 1
                a = Getvcfinfo(each_line2)
                with open(outfile, 'a+') as f3:
                    if flag == 1:
                        print('SVID', 'DatasetID', file = f3)
                        flag = 0
                    datasetid1, datasetid3, region1, region3, gaptype1, gaptype3 = calRegion(t,a,each_line2,outfile, bp)  
               
                        
                    if (len(datasetid1) == 1) or (len(datasetid3) == 1):
                        if region1 != '':
                            if gaptype1[0] not in ['centromere', 'telomere']:
                                gaptype1[0] = 'N'
                            print(a.ID, datasetid1, file = f3)
                            count_vcf += 1
                        elif region3 != '' and  region1 != region3:  
                            if gaptype3[0] not in ['centromere', 'telomere' ]:
                                gaptype3[0] = 'N'                            
                            print(a.ID, datasetid3, file = f3)
                            count_vcf += 1
                            
                    elif (len(datasetid1)) or (len(datasetid3)):
                        if region1 != '':
                            for i in range(len(datasetid1)):
                                if gaptype1[i] not in  ['centromere','telomere']:
                                    gaptype1[i] = 'N' 
                                    gaptype2.append(gaptype1[i])
                                    datasetid2.append(datasetid1[i])
                            print(a.ID, datasetid2, file = f3)
                            gaptype2 = []
                            datasetid2 = []
                        elif region3 != '' and  region1 != region3:
                            for i in range(len(datasetid3)):
                                if gaptype3[i] not in  ['centromere','telomere']:
                                    gaptype3[i] = 'N' 
                                    gaptype4.append(gaptype3[i])
                                    datasetid4.append(datasetid3[i])
                            print(a.ID, datasetid4, file = f3)
                            gaptype4 = []
                            datasetid4 = []
                           
if __name__ == "__main__":
    main() 

    
