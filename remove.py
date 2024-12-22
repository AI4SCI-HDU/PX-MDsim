import PXLink as XL

# 初始化 GromacsSys 对象
topfile = "temp_20_big.top"  # 替换为您的拓扑文件路径
grofile = "temp_20_big.gro"  # 替换为您的坐标文件路径
logfile = "remove_unreacted.log"  # 日志文件
ndxfile = ""  # 如果没有 .ndx 文件，可以留空

# 初始化系统对象
Sys = XL.GromacsSys(topfile=topfile, grofile=grofile, logfile=logfile, ndxfile=ndxfile)

# 移除未反应单体（根据残基名称调整）
unreacted_tma = Sys.remove_residue('TMA')  # 移除未反应的 TMA 单体
unreacted_mpd = Sys.remove_residue('XLN')  # 移除未反应的 MPD 单体

# 输出移除结果
print(f"Removed {unreacted_tma} unreacted TMA residues.")
print(f"Removed {unreacted_mpd} unreacted XLN residues.")

# 导出更新后的系统
Sys.output_contents(new_top="updated.top", new_gro="updated.gro")
