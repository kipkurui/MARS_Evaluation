import os
from os.path import expanduser

# import matplotlib

# matplotlib.use('Agg')
# from Cython.Includes.numpy.__init__ import len
from matplotlib import pyplot as plt

plt.switch_backend('agg')

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# def plot_centrimo(centrimo_in, figure_output):
#     # mot_file = pd.read_table(centrimo_in, index_col=0)
#     centrimo_in.sort(columns="Average", axis=0, ascending=False, inplace=True)
#     cg = sns.clustermap(centrimo_in, method='single', metric="euclidean",
#                         z_score=None, row_cluster=False, col_cluster=True, linewidths=.15, annot=True)
#     test = plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
#     test = plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
#     f = plt.gcf()
#     f.savefig(figure_output, bbox_inches='tight')


def extract_meme(meme_in, motif, meme_out, raw_dict):
    with open(meme_in) as f1:
        with open(meme_out, "a") as out_meme:
            lines = f1.readlines()
            for i, line in enumerate(lines):
                head = line.split()
                if motif in line and motif == head[1]:
                    k = i
                    out_meme.write("\n" + lines[i].strip() + "_" + raw_dict[motif] + "\n\n")
                    if "log-odds" in lines[i + 1]:
                        odds = lines[k + 2].split()
                        for j in range(2, (int(odds[5]) + 3) * 2):
                            out_meme.write(lines[i + j]),
                    elif "letter-probability" in lines[i + 2]:
                        out_meme.write(lines[i + 2])
                        odds = lines[k + 2].split()
                        for j in range(0, (int(odds[5]))):
                            out_meme.write(lines[i + 3 + j])
                    else:
                        print"No data"


def get_dict(raw_file):
    from collections import OrderedDict
    raw_dict = OrderedDict()
    with open(raw_file) as raw_in:
        for line in raw_in:
            raw_dict[line.split()[0]] = line.split()[-1]
    return raw_dict


def get_dict_assess(raw_file, sort_by):
    if sort_by == "MNCP":
        pos = 2
    else:
        pos = 1
    from collections import OrderedDict
    from operator import itemgetter
    # raw_dict = OrderedDict()
    test_d = {}
    with open(raw_file) as raw_in:
        for line in raw_in:
            test_d[line.split()[0]] = line.split()[pos]
    raw_dict = OrderedDict(sorted(test_d.items(), key=itemgetter(1), reverse=True))
    return raw_dict


def extract_scored_meme(meme_in, out_meme, raw_dict):
    print "And finally to extract meme"
    meme_head = """MEME version 4.4\n\nALPHABET= ACGT\n\nstrands: + -\n
    Background letter frequencies (from uniform background):
    A 0.25000 C 0.25000 G 0.25000 T 0.25000\n"""
    with open(out_meme, "w") as meme_out:
        meme_out.write(meme_head)
    for key in raw_dict.iterkeys():
        extract_meme(meme_in, key, out_meme, raw_dict)


def tab2fasta(posneg, fasta, background):
    i = 0
    # print fasta
    with open(posneg) as tab:
        with open(fasta, 'w') as fa:
            with open(background, 'w') as bg:
                for line in tab:
                    details = line.split()
                    if len(details) == 2:
                        pos = 1
                    else:
                        pos = 2
                    if i < 500:
                        fa.write(">" + line.split()[0] + '\n' + line.split()[pos] + "\n")
                    else:
                        bg.write(">" + line.split()[0] + '\n' + line.split()[pos] + "\n")
                    i += 1


def meme2gimme(meme, gimme):
    with open(meme) as motif:
        with open(gimme, 'w') as gmotif:
            for line in motif:
                if line.startswith("MOTIF"):
                    gmotif.write(">" + line.split(" ")[1] + '\n')
                elif line.startswith('letter-probability'):
                    continue
                elif line.startswith('  '):
                    a = line.split()
                    gmotif.write(a[0] + '\t' + a[1] + '\t' + a[2] + '\t' + a[3] + '\n')
                elif line.startswith("\n"):
                    continue
                else:
                    continue


meme_head = '''MEME version 4.4\n
ALPHABET= ACGT\n\nstrands: + -\n
Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 \n'''


def write2meme(meme_out, details):
    outfile = open(meme_out, "a")
    outfile.write(details)


def mkdir_p(path):
    import os
    import errno

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def get_path(method=""):
    home = expanduser("~")
    if home.split('/')[2] == 'kipkurui':
        static_file = '%s/MATOM/static/files/%s' % (BASE_DIR, method)
        meme_path = "%s/meme/bin" % home
    else:
        home = '/home/caleb'
        static_file = '%s/www/MARS/static/files/%s' % (home, method)
        meme_path = "%s/anaconda2/envs/MARS/bin" % home

    if method != "":
        if os.path.isdir(static_file):

            file_list = os.listdir(static_file)

            file_list = map(int, file_list)
            file_list.sort()
        else:
            mkdir_p(static_file)
            file_list = []

        if len(file_list) > 0:
            last_dir = file_list[-1] + 1
        else:
            last_dir = 1

        static_files_path = "%s/%s" % (static_file, str(last_dir))
    else:

        static_files_path = static_file
    return static_files_path, meme_path


meme_path = get_path()[1]


def create_parameter_file(f_name, tf, mode, score=""):
    with open(f_name, 'w') as f_n:
        f_n.write("tf:%s\nmode:%s\nscore:%s\n" % (tf, mode, score))


def get_parameter_dict(f_name):
    dicts = {}
    with open(f_name) as par_file:
        for line in par_file:
            line_det = line.strip().split(":")
            dicts[line_det[0]] = line_det[1]
    return dicts


def combine_meme(meme_in, meme_out):
    found = False
    with open(meme_out, 'a') as meme_out:
        with open(meme_in) as meme:
            for line in meme:
                if len(line.split()) > 1:
                    if line.split()[0] == "MOTIF":
                        found = True
                        meme_out.write("\n\n")
                if found:
                    meme_out.write(line)


def handle_uploaded_bed(tf, length, uploaded_chipseq, results_folder):
    bed_is_correct = is_bed_uploaded(uploaded_chipseq, results_folder)
    if bed_is_correct[0]:
        os.system("%s/MARS_Suite/bed2chipseg.sh %s/%s %s/%s/%s_tmp_chip_%s %s/Data/hg19.fa" %
                  (BASE_DIR, results_folder, uploaded_chipseq, results_folder, tf, tf, length, BASE_DIR))
        return True, "All good"
    else:
        bed_error = bed_is_correct[1]
        return False, bed_error


def get_table(file_in):
    """

    :return:
    """
    gomer_list = []
    with open(file_in) as gomer:
        for line in gomer:
            gomer_list.append(line.split("\t"))
    t_head = gomer_list[0]
    gomer_list = gomer_list[1:]

    return t_head, gomer_list


##############################################################################
# Validation Scripts
###############################################################################


def is_bed_uploaded(bed_in, results_folder):
    import subprocess
    with open('%s/%s' % (results_folder, bed_in.name), 'w') as bed_raw:
        for chunk in bed_in.chunks():
            bed_raw.write(chunk)
    with open('%s/error.txt' % results_folder, 'w') as error:
        subprocess.call('sortBed -i %s/%s >%s/trash' % (results_folder, bed_in.name, results_folder),
                        stderr=error, shell=True)
    if os.stat('%s/error.txt' % results_folder).st_size > 0:
        with open('%s/error.txt' % results_folder, 'r') as output:
            bed_test = output.readlines()
        return False, bed_test
    else:
        return True, "All good"


def is_meme_uploaded(meme_in, results_folder):
    import subprocess
    with open('%s/%s' % (results_folder, meme_in.name), 'w') as meme_raw:
        for chunk in meme_in.chunks():
            meme_raw.write(chunk)

    with open('%s/error.txt' % results_folder, 'w') as error:
        subprocess.call('%s/meme2meme %s/%s >%s/trash' % (meme_path, results_folder, meme_in.name, results_folder),
                        stderr=error, shell=True)

    if os.stat('%s/error.txt' % results_folder).st_size > 0:
        with open('%s/error.txt' % results_folder, 'r') as output:
            meme_test = output.readlines()
        return False, meme_test
    else:
        return True, "All good"


def is_meme_pasted(meme_in, results_folder):
    import subprocess
    with open("%s/pasted.meme" % results_folder, 'w') as meme_raw:
        meme_raw.write(meme_in)
    with open('%s/error.txt' % results_folder, 'w') as error:
        subprocess.call('%s/meme2meme %s/pasted.meme >%s/trash' % (meme_path, results_folder,
                                                                   results_folder), stderr=error, shell=True)
    if os.stat('%s/error.txt' % results_folder).st_size > 0:
        with open('%s/error.txt' % results_folder, 'r') as output:
            meme_test = output.readlines()
        return False, meme_test
    else:
        return True, "All good"


def is_bed(a):
    """
    Quick amd dirty testing if file is bed file

    """
    import types
    number_types = (types.IntType, types.LongType, types.FloatType, types.ComplexType)
    try:
        return isinstance(float(a), number_types)
    except Exception:
        return False
