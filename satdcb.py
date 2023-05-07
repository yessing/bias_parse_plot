'''
deduct sat DCB ref. to OSB files
'''
import math as math

CL = 299792458
SATLIST = [
    'G01','G02','G03','G04','G05','G06','G07','G08','G09','G10',
    'G11','G12','G13','G14','G15','G16','G17','G18','G19','G20',
    'G21','G22','G23','G24','G25','G26','G27'      ,'G29','G30',
    'G31','G32',
    'E01','E02','E03','E04','E05','E06','E07','E08','E09','E10',
    'E11','E12','E13','E14','E15','E16','E17','E18','E19','E20',
    'E21','E22','E23','E24','E25','E26','E27','E28','E29','E30',
    'E31','E32','E33','E34','E35','E36',
    # 'C19','C20','C21','C22','C23','C24','C25','C26','C27','C28',
    # 'C29','C30','C31','C32','C33','C34','C35','C36','C37','C38',
    # 'C39','C40','C41','C42','C43','C44','C45','C46','C47','C48',
    # 'C49','C50',
    # 'R01','R02','R03','R04','R05','R06','R07','R08','R09','R10',
    # 'R11','R12','R13','R14','R15','R16','R17','R18','R19','R20',
    # 'R21','R22','R23','R24','R25','R26','R27','R28','R29','R30',
    # 'R31','R32'
]

chanelDict = {
    'G':['C1C','C2W','C5Q'],
    'E':['C1C','C5Q','C7X'],
    'C':['C2I','C7I','C6I']
}
phchanelDict = {
    'G':['L1C','L2W','L5I'],
    'E':['L1C','L5Q','L7Q'],
    'C':['L2I','L7I','L6I']
}
freqDict = {
    'G':[1e6*1575.42, 1e6*1227.60,1e6*1176.45],
    'E':[1e6*1575.42, 1e6*1176.45,1e6*1207.14],
    'C':[1e6*1561.098,1e6*1207.14,1e6*1268.52]
}
def alpha12(sat:str,pow = 2):
    f1 = freqDict[sat[0]][0]
    f2 = freqDict[sat[0]][1]
    return f1**pow / (f1 * f1 - f2 * f2)    

def beta12(sat:str,pow = 2):
    f1 = freqDict[sat[0]][0]
    f2 = freqDict[sat[0]][1]
    return f2**pow / (f1 * f1 - f2 * f2)
def ionfact(sat:str):
    return 40.3E16 / (freqDict[sat[0]][0]) ** 2

def parselineOSB(line:str,valpos = -2):
    items = line.strip().split()
    epv = int(items[3][9:])     #epoch valid
    return {'sat':items[1],'epv':epv,'chn':items[2],'val':float(items[valpos])} # sat/Chn could be 2，2 


# filepath = "F:\\test_apps\\PPPgamp2e\\GAMP-TEST\\2021.334\\"\
#     +"CAS0MGXRAP_20213340000_01D_01D_OSB.BIA"
filepath = "F:\\data\\bia\\cnt22415.bia"

def getsdcb():
    fo = open(filepath)
    allrec = {}
    alldcbrec = {}

    for line in fo:
        if line[1:4] == "OSB":
            rec = parselineOSB(line)
            if not rec['sat'] in SATLIST:
                continue
            if not rec['sat'] in allrec:
                allrec[rec['sat']] = {}
            allrec[rec['sat']][rec['chn']] = rec['val']

    for isat in allrec:
        allrec[isat]['dcb'] = (allrec[isat][chanelDict[isat[0]][0]] - \
            allrec[isat][chanelDict[isat[0]][1]]) * CL * 1e-9 * beta12(isat)/\
            ionfact(isat)
        alldcbrec[isat] = allrec[isat]['dcb']   # unit TECU
        # print(isat + "  "+ str(allrec[isat]['dcb']))

    # print(len(allrec))
    return alldcbrec

def getsdnl():
    fo = open(filepath)
    allrec = {}
    alldnlrec = {}
    alleprec = {}
    allepdnlrec = {}
    ep_last = -1
    ep_now = 0

    for line in fo:
        if line[1:4] == "OSB":
            rec = parselineOSB(line,-1)     # varpos=-1 if file contains no std
            if not rec['epv'] in alleprec:
                alleprec[rec['epv']] = {}

            if not rec['sat'] in SATLIST:
                continue
            if not rec['sat'] in alleprec[rec['epv']]:
                alleprec[rec['epv']][rec['sat']] = {}
            alleprec[rec['epv']][rec['sat']][rec['chn']] = rec['val']

    for iep in alleprec:
        allrec = alleprec[iep]
        alldnlrec = {}
        for isat in allrec:
            try:
                sysf = isat[0]
                _f1 = freqDict[sysf][0]
                _f2 = freqDict[sysf][1]     # unit Hz
                _b1 = allrec[isat][chanelDict[sysf][0]] * 1e-9
                _b2 = allrec[isat][chanelDict[sysf][1]] * 1e-9    # unit second
                _B1 = allrec[isat][phchanelDict[sysf][0]] * 1e-9 * _f1
                _B2 = allrec[isat][phchanelDict[sysf][1]] * 1e-9 * _f2    #  unit cyc

                allrec[isat]['dnl'] = (_f1 * _B1 - _f2 * _B2) / (_f1 - _f2) -\
                     (_f1 * _f1 * _b1 - _f2 * _f2 * _b2) / (_f1 - _f2)
                alldnlrec[isat] = allrec[isat]['dnl']
                alldnlrec[isat] = _b2*CL*100
                # alldnlrec[isat] = allrec[isat]['dnl'] - round(allrec[isat]['dnl'])   # only record fractional part
            except:
                print(f"key error check input file epoch = {iep} sat = {isat}\n")
                alldnlrec[isat] = math.nan
            # print(isat + "  "+ str(allrec[isat]['dcb']))
        alleprec[iep].clear()
        allepdnlrec[iep] = alldnlrec
    
    del alleprec
    # print(len(allrec))
    return allepdnlrec    

from matplotlib import pyplot as plt
if __name__ == "__main__":
    res = getsdnl()
    # rkey = list(res.keys())[-1]
    # res = res[rkey]
    # dfrac = {}
    subsatli = {'G02':[],'G05':[],'G10':[],'G12':[],'G15':[],
        'G18':[],'G23':[],'G24':[],'G25':[],'G29':[]}
    for iep in res:
        for isat in subsatli:
            try:
                subsatli[isat].append(res[iep][isat])
            except:
                subsatli[isat].append(math.nan)
                print(f"plot data missing at {iep} {isat}")
    for isat in subsatli:
        plt.plot(subsatli[isat],label=isat)

    plt.title("C2W OSB",loc="left")
    plt.ylabel("bias/cm")
    plt.xlabel("epoch/$\\times 30s$")
    plt.legend()
    plt.grid()
    plt.show()
    print("done!\n")
