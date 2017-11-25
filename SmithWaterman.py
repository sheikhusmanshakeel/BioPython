#!/usr/bin/python
from math import *

strG = str(raw_input("Genome sequence: "))
strR = str(raw_input("Read sequence: "))
strG = strG.upper()
strR = strR.upper()
print "Your genome ref input : " + strG
print "Your read ref input   : " + strR

matrix = []
path = []
row = len(strR)
col = len(strG)
strR = "^" + strR
strG = "^" + strG
for i in range(row + 1):
    matrix.append([0] * (col + 1))
    path.append(["N"] * (col + 1))


def print_matrix(matrix):
    print '\t' + ('\t'.join(map(str, list(strG))))
    i = 0
    for line in matrix:
        print strR[i] + "\t" + ('\t'.join(map(str, line)))
        i += 1


# print_matrix(matrix)
indelValue = -1
matchValue = 2
for i in range(1, row + 1):
    for j in range(1, col + 1):
        # penalty map
        from_left = matrix[i][j - 1] + indelValue
        from_top = matrix[i - 1][j] + indelValue
        if strR[i] == strG[j]:
            from_diag = matrix[i - 1][j - 1] + matchValue
        else:
            from_diag = matrix[i - 1][j - 1] + indelValue

        matrix[i][j] = max(from_left, from_top, from_diag)
        # path map
        if matrix[i][j] == from_left:
            path[i][j] = "-"
        elif matrix[i][j] == from_top:
            path[i][j] = "|"
        elif matrix[i][j] == from_diag:
            path[i][j] = "M"
        else:
            pass

        if matrix[i][j] < 0:
            matrix[i][j] = 0
pass
print_matrix(matrix)
print
print_matrix(path)

# trace back

iRow = len(matrix) - 1
jCol = len(matrix[0]) - 1

# def report_trace_point(x,y,strX,strY,pathMap,reportX,reportY):
# 	if pathMap[x][y] == "M":
# 		reportX.append(strX[x-1])
# 		reportY.append(strY[y-1])
# 		report_trace_point([x-1,y-1,strX,strY,pathMap,reportX,reportY])
# 		# return report_trace_point( x-1, y-1, strX,strY,pathMap)
# 	elif pathMap[x][y]=="-":
# 		reportX.append(strX[x])
# 		reportY.append(strY[y-1])
# 		report_trace_point([x,y-1,strX,strY,pathMap,reportX,reportY])
# 		# return report_trace_point( x,y-1 , strX,strY,pathMap)
# 	elif pathMap[x][y]=="|":
# 		reportX.append(strX[x-1])
# 		reportY.append(strY[y])
# 		report_trace_point([x-1,y,strX,strY,pathMap,reportX,reportY])
# 		# return report_trace_point(x-1 ,y , strX,strY,pathMap)
# 	else:
# 		reportX.append(strX[x])
# 		reportY.append(strY[y])
# 		report_trace_point([x,y,strX,strY,pathMap,reportX,reportY])
# 		# return report_trace_point( x,y , strX,strY,pathMap)

# print report_trace_point(4,3,strR,strG,path)
while iRow >= 0:
    maxPnltyVaue = max(matrix[iRow])
    while jCol >= 0:
        if matrix[iRow][jCol] == maxPnltyVaue:
            ipair = iRow
            jpair = jCol
            reportR = []
            reportG = []
            while (1):
                if ipair == 0 and jpair == 0:
                    break
                # else:
                if path[ipair][jpair] == "M":
                    reportR.append(strR[ipair])
                    reportG.append(strG[jpair])
                    ipair -= 1
                    jpair -= 1
                elif path[ipair][jpair] == "-":
                    # reportR.append(strR[ipair])
                    reportR.append("-")
                    reportG.append(strG[jpair])
                    # ipair -= 1
                    jpair -= 1
                elif path[ipair][jpair] == "|":
                    reportR.append(strR[ipair])
                    reportG.append("-")
                    ipair -= 1
                # jpair -= 1
                elif path[ipair][jpair] == "N":
                    if ipair > 0:
                        reportR.append(strR[ipair])
                        ipair -= 1

                    if jpair > 0:
                        reportG.append(strG[jpair])
                        jpair -= 1

            print "Your alignment would be followed: "
            print "R:" + "".join(reversed(reportR))
            print "G:" + "".join(reversed(reportG))
            print
        jCol -= 1
    iRow -= 1
