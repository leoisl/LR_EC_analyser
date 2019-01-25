def fromStringToRoundedValues(s):
	return [round(float(x), 2) for x in s.split()]

def calculate(deletions, homDel):
	for d, h in zip(deletions, homDel):
		print(round(h/d*100, 1))


deletions=fromStringToRoundedValues("7.41 0.85 0.64 0.17 2.15 1.73 0.39 2.64 1.94 0.09 0.3 1.51 0.18  4.82 3.82 3.27 2.51 4.08 5.06 5.17 ")
homDeletions=fromStringToRoundedValues("2.96 0.28 0.19 0.03 0.77 0.63 0.19 0.62 0.42 0.02 0.1 0.46 0.04  2.46 2.14 2.05 1.82 2.05 2.26 2.27 ")
calculate(deletions, homDeletions)


insertions=fromStringToRoundedValues("1.2 0.21 0.14 0.05 0.32 0.24 0.07 0.47 0.4 0.08 0.19 0.22 0.03 0.28 0.28 0.19 0.04 0.06 0.09 0.05 ")
homInsertions=fromStringToRoundedValues("0.38 0.05 0.03 0.01 0.09 0.07 0.02 0.11 0.09 0.01 0.02 0.06 0.01 0.08 0.06 0.03 0.01 0.01 0.02 0.01 ")
calculate(insertions, homInsertions)