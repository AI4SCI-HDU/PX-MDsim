# 获取整个系统中所有苯环的法向量，据此求有序系数和密度分布。
# 自制程序，直接读取.gro文件。
# 有序系数为S(r) = <(3cos^2(θ(r)) - 1) / 2>
# 密度分布为G(r) = ρ(r) / ρ(bulk)
# 程序读取xtc文件，由subprocess调用gmx trjconv -dump来生成多个
# 对应同不同时刻的gro文件，分别计算结果，并输出平均结果。
# ###### Changelog #######
# 22/12/08 更新：记录并输出距离过近的苯环对；修正了导致计数两次的bug
# 22/12/09 更新：修正了对“XY中有一方向跨过PBC，另一方向位于正中间”的苯环的处理
# 22/12/10 更新：再次修正导致计数两次的bug
# 22/12/12 更新：修改了r=rmax时的处理

import numpy as np
import regex as re
import math
import subprocess as sp
import os


# 定义储存数据的class
class ring_residue:
    """存储含苯环残基相关信息的class
    包含如下变量：
    id: 整数，记录res编号
    name: 字符串，记录残基名
    ar_carbon: 6*3的numpy array，记录6个芳香碳的座标；
    center: numpy array，记录苯环中心的座标；
    normal: numpy array，记录苯环法向量的矢量；
    """

    def __init__(self, resid, resname, carbon_coord):
        self.id = resid
        self.name = resname
        self.ar_carbon = np.array(carbon_coord)
        # 若芳香碳原子数不为6，报错
        if len(self.ar_carbon) != 6:
            now_ar_carbons = self.ar_carbon[0:]  # 获取芳香碳原子
            print(f"carbons for residue {self.id}: {now_ar_carbons}")
            raise ValueError(
                "There are supposed to be 6 aromatic carbon atoms.")
        # 计算苯环中心
        self.center = np.mean((self.ar_carbon), axis=0)

    def fix_pbc_skip(self, box_dim):
        """部分苯环跳过了周期性边界截，导致其原子位置相差巨大。
        需要设法将原子移至同一侧，才能正常计算中心点和法向量。
        """
        # 若苯环被周期性边界截断，则苯环内会有XY方向上距离巨大的原子。记录这些残基。
        atom_dist_XY = self.ar_carbon.max(axis=0) - self.ar_carbon.min(axis=0)
        pbc_jump_XY = [False, False]
        if (atom_dist_XY[0] * 2 > box_dim[0]):
            pbc_jump_XY[0] = True
        if (atom_dist_XY[1] * 2 > box_dim[1]):
            pbc_jump_XY[1] = True
        pbc_jump = pbc_jump_XY[0] or pbc_jump_XY[1]
        # 为了防止“XY中有一方向跨过PBC，另一方向位于中线”时出现问题，需对X、Y方向单独分别处理
        for dim in range(2):
            if pbc_jump_XY[dim]:
                # 首先将各原子都放回box内部。
                # 由于我们的体系准备，残基不可能在Z方向上超出周期性边界，故暂且不考虑。
                # 我原本准备采用取余的方法将座标放回周期性边界：
                # self.ar_carbon[:, 0] = self.ar_carbon[:, 0] % box_dim[0]
                # self.ar_carbon[:, 1] = self.ar_carbon[:, 1] % box_dim[1]
                # 但考虑到可能出现的浮点数误差，还是改用只需加减法的方法。
                # 同时，记录是否有原子在这一步发生移动。
                pbc_shift = False
                for i in range(6):
                    while (True):
                        if self.ar_carbon[i][dim] < 0:
                            self.ar_carbon[i][dim] += box_dim[dim]
                            pbc_shift = True
                        elif self.ar_carbon[i][dim] > box_dim[dim]:
                            self.ar_carbon[i][dim] -= box_dim[dim]
                            pbc_shift = True
                        else:
                            break
                # 检查每个原子各位于周期性边界的哪一边。
                # [0]记录接近0一侧的原子，[1]记录接近最大值一侧的原子
                sep = [[], []]
                for i in range(6):
                    if self.ar_carbon[i][dim] <= box_dim[dim] / 2:
                        sep[0].append(i)
                    else:
                        sep[1].append(i)
                # 移动原子数量较少的一侧的原子
                if len(sep[0]) <= len(sep[1]):
                    to_move = sep[0]
                    # 将接近0一侧移至接近最大值一侧，x' = x + box
                    move_by = 1
                else:
                    to_move = sep[1]
                    # 将接近最大值一侧移至接近0一侧，x' = x - box
                    move_by = -1
                if to_move:
                    self.ar_carbon[to_move, dim] += box_dim[dim] * move_by
                # 存在一种状况，仅需把部分超出周期性边界的原子放回box内部，就能完成修正。
                # 这时会出现pbc_jump=True且to_move为空的状况。
                # 检查修正失败的状况，即“没有超过周期性修正，没有需要移动的原子，但pbc_jump=True”
                if (not pbc_shift) and (not to_move):
                    raise ValueError(
                        f"Falied in PBC fixing of benzene ring in residue {self.id}."
                    )
        # 重新计算中心位置
        if pbc_jump:
            self.center = np.mean((self.ar_carbon), axis=0)
        return pbc_jump

    def calc_norm(self):
        """计算苯环法向量，存入self.normal变量中。
        我们需要先计算苯环中心位置，判断它是否在我们关注的位置范围内，若是，才有可能用到其法向量。
        出于优化目的，将这部分从__init__中独立出来，按需调用。

        """
        # 计算苯环法向量
        vector_a = self.center - self.ar_carbon[0]
        vector_ab_parallel = True
        # 寻找一对不平行的中心-芳香碳矢量。
        for i in range(1, 6):
            vector_b = self.center - self.ar_carbon[i]
            norm_vector = np.cross(vector_a, vector_b)
            if np.linalg.norm(norm_vector) > 0.000001:
                vector_ab_parallel = False
                break
        if vector_ab_parallel:
            raise ValueError(
                "Something's wrong, we can't find two non-parallel center - carbon vectors."
            )
        self.normal = norm_vector


def pbc_dist(coord_A, coord_B, box_dim):
    """计算A，B两点在周期性边界box_dim中的距离

    Args:
        coord_A (numpy.array): 点A座标
        coord_B (numpy.array): 点B座标
        box_dim (numpy.array): 周期性边界大小
    """
    distance = np.abs(coord_A - coord_B)
    distance_pbc = np.where(distance > 0.5 * box_dim, box_dim - distance,
                            distance)
    return np.sqrt((distance_pbc**2).sum())


def cos_angle(vector_A, vector_B):
    """A，B两个向量夹角为θ，计算cosθ

    Args:
        vector_A (numpy.array): 向量A
        vector_B (numpy.array): 向量B
    """
    cos_theta = np.dot(
        vector_A,
        vector_B) / np.linalg.norm(vector_A) / np.linalg.norm(vector_B)
    return cos_theta


# 将计算部分划分到独立的函数中
def rings_param(gro):
    """读取输入的gro文件，寻找其中的苯环，计算径向分布函数G(r)和有序参数S(r)
    寻找苯环需要的残基名、原子名及计算所用的Z轴范围由程序开头由程序开头给出的全局变量指定。

    Args:
        gro (string): 输入gro文件的路径

    Raises:
        ValueError: gro文件格式不对时，报错
        ValueError: 某一距离上苯环对数为0但有序参数S(r)不为0时，报错

    Returns:
        list: 三个不同的list，分别对应不同距离上的苯环对数，径向分布函数G(r)和有序参数S(r)
    """
    # 判断残基编号是否发生变化，同时也作为读取第一个残基时的处理flag
    current_resid = 0
    current_resname = ''
    # 储存当前残基的各个芳香碳位置座标
    current_coord = []
    # 记录残基数
    num_res = 0
    num_ring_sel = 0
    num_ring_ref = 0
    # 记录用的list。这一list将包含每个“含苯环残基”对应的ring_residue类。
    residues = []
    # 记录哪些苯环位于ref和sel范围内，内容为对应residue在residues list内的序号。
    ref_ix = []
    sel_ix = []
    # 记录跃过周期性边界的残基。这一list的内容为这些残基对应的residues list序号，用于输出、检查结果。
    skipped_residues = []
    # 记录距离在指定watch范围内的苯环对。
    marked_ring = []
    # 计数用的array
    pairnum_ring = np.zeros(n_bin)
    order_ring = np.zeros(n_bin)
    with open(gro, mode='r') as f:
        # 跳过文件标题。不需要读取时间和总原子数。
        f.readline()
        f.readline()
        # 读取关于原子数据的第一行
        # 逐行读取gro文件，直到到达文件末尾。
        while (True):
            line = f.readline().split()
            # 文件末行为box尺寸。对于正常的.gro文件，仅有在此时，
            # 才有len(line) == 3。此时，读取box尺寸并中止循环。
            if len(line) == 3:
                box_dim = np.array([float(i) for i in line])
                break
            # 判断.gro文件是否含有速度
            if current_resid == 0:
                if len(line) == 6:
                    # 当.gro文件不包括速度时，一行数据.split()后有6项，倒数三项为座标
                    coord_ix = [-3, -2, -1]
                elif len(line) == 5:
                    # 当.gro文件不包括速度时，一行数据.split()后有5项，倒数三项为座标
                    coord_ix = [-3, -2, -1]
                elif len(line) == 9:
                    # 当.gro文件包括速度时，一行数据.split()后有9项，倒数4~6项为座标
                    coord_ix = [-6, -5, -4]
                elif len(line) == 8:
                    # 当.gro文件包括速度时，一行数据.split()后有9项，倒数4~6项为座标
                    coord_ix = [-6, -5, -4]
            # 读取残基号和残基名，两者均在line[0]中。残基号为数字，残基名为字母或字母+数字
            res_split = re.findall(r'[A-Za-z]+|\d+', line[0])
            resid = int(res_split[0])
            if len(res_split) == 3:
                resname = res_split[1] + res_split[2]
            else:
                resname = res_split[1]
            # 判断是否是新残基。如是，则修改“当前残基”的编号
            if resid != current_resid:
                # 且如当前残基并非第一个残基(current_resid不为0)，且上一个残基是含苯环残基（current_coord不为空）
                # 则上一残基的数据计算完毕，将其记录。
                if current_resid and current_coord:
                    current_residue = ring_residue(current_resid,
                                                   current_resname,
                                                   current_coord)
                    residues.append(current_residue)
                    if ref_zmin <= current_residue.center[2] <= ref_zmax:
                        # 如果这一苯环处在{zmin, zmax}中间，则将其计入list residues
                        ref_ix.append(num_res)  # 当录入第一个环时，num_res=1，故需要减1
                        num_ring_ref = num_ring_ref + 1
                    elif sel_zmin <= current_residue.center[2] <= sel_zmax:
                        # sel_ix记录“只”在sel范围内而不在ref范围内的苯环
                        sel_ix.append(num_res)
                        num_ring_sel = num_ring_sel + 1
                    num_res = num_res + 1
                current_coord = []
                current_resid = resid
                current_resname = resname
            # 判断该残基是否是含有苯环的残基
            if any([(resn in resname) for resn in resname_ring]):
                # 如是，则开始留意原子名
                # 当原子编号不多于四位数时，原子名是line[1]
                if len(line) == 6 or len(line) == 9:
                    atom_name = line[1]
                # 否则，原子名要从line[1]中提取。这里取前三位。实际上这个方法对SOL不适用
                # (水中氧的原子名为OW)，但我们本来就不关心水，且resname == SOL时，这一段不会执行。
                elif len(line) == 5 or len(line) == 8:
                    atom_name = line[1][:-5]
                else:
                    # 根据.gro文件是否记录速度及原子编号是否为5位，原子行可能有5，6，8，9段。
                    # 如果len(line)为其他值，则输入文件有问题，报错。
                    print("Error at line 111\n")
                    raise ValueError(
                        "Something's wrong with the input .gro file.")
                # 如果这个原子属于芳香碳，则记录其座标，用于计算。
                if atom_name in ar_carbon:
                    atom_coord = np.array([float(line[i]) for i in coord_ix])
                    current_coord.append(atom_coord)

    num_ring_all = num_ring_sel + num_ring_ref
    # 读取gro文件的最后一步读取了box的大小。用这一信息检查哪些苯环跃过了周期性边界，
    # 对其进行修正，然后计算法向量。
    for i in residues:
        skipped = i.fix_pbc_skip(box_dim)
        if skipped:
            skipped_residues.append(i.id)
        # 记录哪些残基经过了这种修正。
        i.calc_norm()

    # 循环每一对苯环
    for i in range(num_ring_ref):
        # 首先循环ref-ref苯环对。此时需要注意防止重复计数。
        for j in range(i + 1, num_ring_ref):
            ri = residues[ref_ix[i]]
            rj = residues[ref_ix[j]]
            # 计算中心距离
            distance_ij = pbc_dist(ri.center, rj.center, box_dim)
            # 如果这两个苯环距离过近，则对其记录。
            if distance_ij < dist_too_low:
                pairs_too_close.append([gro, ri.id, rj.id, distance_ij])
            # 如果这两个苯环间距为指定的watch范围，对其记录。
            if watch_dist and (watch_dist[0] < distance_ij < watch_dist[1]):
                marked_ring.append([ri.id, rj.id, distance_ij])
            bin_ij = math.floor(distance_ij / dr + 0.5)
            # 判断这一距离落在哪个bin内。如果落在距离外，跳过这一对苯环。
            # NOTE: 实际上由于bin_ij的计算方法，n_bin个bin对应的r范围分别是
            # [0, 0.5dr], [0.5dr, 1.5dr], ... [(n_bin - 1.5)dr, (n_bin - 0.5)dr]
            # 实际上最后结果不会反映[(n_bin - 0.5)dr, n_bin * dr]范围内的对。
            if bin_ij >= n_bin:
                continue
            # 计数。
            pairnum_ring[bin_ij] = pairnum_ring[bin_ij] + 1
            # 计算(3cos^2(θ(r)) - 1) / 2
            # 先计算cos(θ(r))
            cos_theta = cos_angle(ri.normal, rj.normal)
            # TODO:这里找到的是ij和ji两对苯环，那么计算结果是否需要相应地乘以2？目前按乘以2算。
            order_ij = (3 * (cos_theta**2) - 1) / 2
            order_ring[bin_ij] = order_ring[bin_ij] + order_ij
        # 然后循环ref-sel苯环对。不会出现重复计数。
        for j in range(num_ring_sel):
            ri = residues[ref_ix[i]]
            rj = residues[sel_ix[j]]
            # 计算中心距离
            distance_ij = pbc_dist(ri.center, rj.center, box_dim)
            # 判断这一距离落在哪个bin内
            if distance_ij > rmax:
                continue
            # 如果这两个苯环距离过近，则对其记录。
            if distance_ij < dist_too_low:
                pairs_too_close.append([gro, ri.id, rj.id, distance_ij])
            bin_ij = math.floor(distance_ij / dr + 0.5)
            if bin_ij == n_bin:
                bin_ij = bin_ij - 1
            # 计数。
            pairnum_ring[bin_ij] = pairnum_ring[bin_ij] + 1
            # 计算(3cos^2(θ(r)) - 1) / 2
            # 先计算cos(θ(r))
            cos_theta = cos_angle(ri.normal, rj.normal)
            order_ij = (3 * (cos_theta**2) - 1) / 2
            order_ring[bin_ij] = order_ring[bin_ij] + order_ij
    # 求体系平均值
    # 待计算区域的体积
    bulk_vol = box_dim[0] * box_dim[1] * (sel_zmax - sel_zmin)
    # 待计算区域内苯环的数量密度，即G(r)公式中的ρ(bulk)
    bulk_ring_density = num_ring_all / bulk_vol
    # 各个半径{r, r+dr}球层的体积，近似为4*pi*r^2*dr
    vol_dr = np.array(
        [4 * math.pi * (dr * (i + 0.5))**2 * dr for i in range(n_bin)])
    # 求G(r)、S(r)
    ring_radial_distribution = pairnum_ring / bulk_ring_density / vol_dr / num_ring_ref
    ring_order_parameter = []
    for i in range(n_bin):
        if pairnum_ring[i] == 0:
            if order_ring[i] != 0:
                raise ValueError(
                    "Somehow we calculated ring order from no ring pairs.")
            else:
                ring_order_parameter.append(0)
        else:
            ring_order_parameter.append(order_ring[i] / pairnum_ring[i])
    return pairnum_ring, ring_radial_distribution, ring_order_parameter, marked_ring


def extract_frame(time):
    """使用subprocess调用gmx trjconv -dump，从xtc文件中提取出指定时刻的一帧。
    gmx trjconv所需的其他参数均为程序开头给出的全局变量。

    Args:
        time (_type_): _description_
    """
    out_gro = frame_file + str(time) + ".gro"
    cmd = f"echo 0 | gmx trjconv -f {xtc} -s {tpr} -o {out_gro} -dump {str(time)}"
    sp.check_call(cmd, shell=True)
    return out_gro


# 指定输入文件
# xtc = 'DPC70/sd70p30f_ext2_pd3part.xtc'
# xtc = 'Mem70/d70_pd2_50ns.xtc'
xtc = 'tm60.xtc'
# tpr文件，仅在gmx trjconv -dump时用到
# tpr = 'DPC70/sd70p30f_ext2_pd3.tpr'
# tpr = 'Mem70/d70_pd2.tpr'
tpr = 'tm60.tpr'
# gmx trjconv -dump输出gro文件的文件名头
# 实际文件名为frame_file + str(time) + ".gro"
# frame_file = 'DPC70/sd70p30f_pd_'
# frame_file = 'Mem70/d70_pd2_'
frame_file = 'sd30_pd_'
# 指定输出文件
# out = 'DPC70/Ring_parameters_sd70_400to500_r1b50_2.txt'
# out = 'Mem70/Ring_parameters_mem70_150to200_r1b50.txt'
out = 'Ring_parameters_sd30_pd5_z1r2b100.txt'
# 程序可以寻找距离在一定范围内的苯环对，并将其输出。
# 为了便于进一步分析，输出时每帧输出一个文件，仅输出苯环残基编号对，便于直接由程序读取。
# 文件名为以下字符串指定的开头+帧时间+.txt
# out_pairs = "DPC70/sd70_pairs_t"
# out_pairs = "Mem70/mem70_pairs_t"
out_pairs = "sd30_pairs_t"
# 指定开始时间、结束时间、时间点数，单位ps
# b_time = 400000
# e_time = 500000
# ntime = 21
b_time = 1500000
e_time = 1550000
ntime = 21
time_points = np.linspace(b_time, e_time, ntime)
dt = time_points[1] - time_points[0]
# 当两个苯环之间距离小于dist_too_low时，将其记录，最终输出。单位nm
dist_too_low = 0.3
pairs_too_close = []
# 监视距离在这一范围内的苯环对，输出其残基编号
watch_dist = []  # watch_dist应为包含两个值的list，分别为监视范围的最短和最长距离
marked_traj = []
# 是否输出单帧结果
output_single_frame = True
# 是否保留dump的gro文件
keep_dump_gro = False
# 补充说明
discription = "DPC 30 solvated membrane, 1500-1550 ns, r=2, L=1"

# 各种残基名。寻找所有包含苯环的残基，即包括TMA、MPD的残基
resname_ring = ["TMA", "TMAm", "TMAd", "TMAt", "XLN", "XLNm", "XLNd"]  # TMA和MPD残基
# 残基中的6个芳香碳原子名
ar_carbon = ["C7", "C8", "C9", "C10", "C0A", "C12"]  # TMA芳香碳原子名
ar_carbon += ["C1", "C2", "C3", "C4", "C5", "C6"]  # MPD芳香碳原子名

# 有序系数和密度分布的最大值和bin数
rmax = 5
n_bin = 500
dr = rmax / n_bin

# Z轴范围的最大最小值，单位nm
# NOTE:因为我们需要算ρ(bulk)，所以必须指定一个Z轴范围，将膜主体部分框在这一范围内，
# 而不能简单地使用整张膜甚至整个box
# 像之前的RDF程序一样，分为ref（中心）和sel（周围）两个范围，不统计sel-sel中的苯环对，
# 确保所有volume都是完整的球层。
# 假定这一范围在所有时刻都不变。
ref_zmin = 1.0
ref_zmax = 6.0
sel_zmin = ref_zmin - rmax
sel_zmax = ref_zmax + rmax

# 各个bin的中点
bin_midpoint = np.array([dr * (i + 0.5) for i in range(n_bin)])

# 记录各个gro文件的计算结果
pairnum_traj = []
distribution_traj = []
order_traj = []

# 对各个文件分别进行计算
for t in time_points:
    f = extract_frame(t)
    pairnum_frame, distribution_frame, order_frame, marked_frame = rings_param(
        f)
    pairnum_traj.append(pairnum_frame)
    distribution_traj.append(distribution_frame)
    order_traj.append(order_frame)
    marked_traj.append(marked_frame)
    # 视需要删除gro文件
    if not keep_dump_gro:
        os.remove(f)

# 求所有帧的平均
pairnum_mean = np.mean(np.array(pairnum_traj), axis=0)
distribution_mean = np.mean(np.array(distribution_traj), axis=0)
order_mean = np.mean(np.array(order_traj), axis=0)
# 输出结果
with open(out, mode="w") as f:
    f.write("Benzene ring radial distribution and order parameter.\n")
    f.write(f"Calculate using coordinate {xtc}\n")
    f.write(f"{discription}\n")
    f.write(
        f"Beginning time {b_time}, end time {e_time}, time interval {dt}\n")
    f.write(f"Reference Z range {ref_zmin} to {ref_zmax}\n")
    f.write(f"Selection Z range {sel_zmin} to {sel_zmax}.\n")
    f.write(f"Average of all {ntime} time points:\n")
    f.write(
        "Bin midpoint\tPair numbers\tRadial distribution\tOrder parameter:\n")
    for i in range(n_bin):
        f.write(
            f"{(bin_midpoint[i]):.3f}\t{pairnum_mean[i]:.2f}\t{distribution_mean[i]:.3f}\t{order_mean[i]}\n"
        )
    if pairs_too_close:
        f.write("Warning: Some pairs are too close.\n")
        f.write("They are listed below: file, index_i, index_j, distance\n")
        for i in pairs_too_close:
            for j in i:
                f.write(str(j) + "\t")
            f.write("\n")
    if output_single_frame:
        for i in range(ntime):
            f.write(f"\nResults calculated at time {time_points[i]}:\n")
            f.write(
                "Bin midpoint\tPair numbers\tRadial distribution\tOrder parameter:\n"
            )
            for j in range(n_bin):
                f.write(f"{(bin_midpoint[j]):.3f}\t{pairnum_traj[i][j]}\t")
                f.write(f"{distribution_traj[i][j]:.3f}\t{order_traj[i][j]}\n")
    if watch_dist:
        f.write(
            f"Found these ring pairs in range {watch_dist[0]:.3f} - {watch_dist[1]:.3f}\n"
        )
        f.write("These pairs are saved in given files.\n")
if watch_dist:
    for i in range(ntime):
        fn = out_pairs + str(time_points[i]) + ".txt"
        with open(fn, 'w') as f:
            for j in marked_traj[i]:
                print(j, file=f)
