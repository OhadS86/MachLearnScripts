import itertools


ind = [];
data = [];
dep = 3; # dep = order + 1
gInfo = []; # Group Information: the detailed combinations of indices
L = 7;   # length of a single sequence

# Generate the variable: gInfo
for i in xrange (dep):
	gInfo.extend (list (itertools.combinations (xrange (L), i + 1)));

fout = open ('gInfo.txt', 'w');
for tp in gInfo:
	fout.write (' '.join (map (str, tp)) + '\n' );
fout.close ();


# Generate the variable: ind
ind.append (0);
for tp in gInfo:
	ind.append (ind[-1] + 4**(len (tp)));
fout = open ('ind.txt', 'w');
fout.write ( ' '.join (map (str, ind)) + '\n');
fout.close ();

# Generate the variable: data, lab
val = {'a':0, 'c':1, 'g': 2, 't':3};
def process (s, tp):
	m = 0;
	for ch in tp:
		m = m * 4 + val[s[ch]];
	ret = [0]*(4**len(tp));
	ret[m] = 1;
	return ret;

fin = open ('testfile5_hs', 'r');
fout = open ('test_data.txt', 'w');
fout_lab = open ('test_lab.txt', 'w');
for line in fin:
	line = line.strip ().lower ();
	if line[0] == '>':
		if line[2] == '5':
			fout_lab.write ('1\n');
		else:
			fout_lab.write ('-1\n');
		continue
	tmpLst = [];
	for tp in gInfo:
		tmpLst.extend (process (line, tp));
	# Replace the following code after installing SciPy and NumPy
	fout.write (' '.join (map (str, tmpLst)) + '\n')

fout.close ();
fin.close ();




