# 文件：run_analysis.py
# 导入我们编写的扩展模块
from gget_pdb import GGETPDB
# 1. 创建分析工具实例
analyzer = GGETPDB()

# 2. 示例：从基因名查找结构
gene_name = "INS"  # 可以替换为 TP53, INS, ACE2 等

structures = analyzer.gene_to_structures(gene_name)


if structures:
    print(f"找到结构: {structures}")

    # 3. 示例：获取第一个结构的详细信息
    first_pdb = structures[0]
    info = analyzer.fetch_pdb_info(first_pdb)
    if info:
        print(f"\n结构 {first_pdb} 的信息:")
        print(f"  标题: {info.get('title')}")
        print(f"  分辨率: {info.get('resolution')}")

    # 4. 示例：生成并打印分析报告
    print("\n" + "=" * 50)
    report = analyzer.generate_report(gene_name=gene_name)
    print(report)

    # 5. 注意：3D可视化需要在Jupyter Notebook中运行
    print("\n提示：如需交互式3D可视化，请将以下代码复制到Jupyter Notebook中运行：")
    print(f"  from gget_pdb import gget_pdb")
    print(f"  viewer = gget_pdb.view_3d('{first_pdb}')")
    print(f"  viewer.show()")

else:
    print(f"未找到基因 '{gene_name}' 的相关结构。")