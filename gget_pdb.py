#文件:gget_pdb.py
import requests
import json
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import ProtParam
import py3Dmol
import pandas as pd
import numpy as np
from IPython.display import display, Markdown, HTML
import warnings

warnings.filterwarnings('ignore')


class GGETPDB:
    """gget的PDB结构分析扩展"""

    def __init__(self):
        self.rcsb_base = "https://data.rcsb.org/rest/v1"
        self.uniprot_api = "https://rest.uniprot.org/uniprotkb"

    # ==================== 1. 智能映射 ====================
    def gene_to_structures(self, gene_name, species="human", max_structures=5):
        """将基因名映射到相关PDB结构"""
        print(f"🔍 正在查询基因 '{gene_name}' 的蛋白结构...")

        # 使用gget获取基因信息
        try:
            import gget
            search_result = gget.search(gene_name, species=species)
            # 正确判断DataFrame是否为空，并提取第一个基因的ID
            if search_result.empty:  # 使用 .empty 属性判断
                return []

            gene_id = search_result.iloc[0]['ensembl_id']
            info_df = gget.info([gene_id])  # 返回的是一个DataFrame
            # 从DataFrame中提取‘uniprot_id’列，如果没有该列则为None
            if not info_df.empty and 'uniprot_id' in info_df.columns:
                uniprot_id = info_df.iloc[0]['uniprot_id']
                # 处理可能存在的多个ID（比如用分号隔开的情况）
                if pd.notna(uniprot_id):
                    # 取第一个ID（如果需要所有ID，可以保留列表）
                    uniprot_id = str(uniprot_id).split(';')[0].strip()
                else:
                    uniprot_id = None
            else:
                uniprot_id = None

            if not uniprot_id:
                # 备用方案：直接通过UniProt API搜索
                params = {"query": f"gene:{gene_name} AND organism:{species}", "format": "json"}
                response = requests.get(self.uniprot_api, params=params).json()
                if response.get("results"):
                    uniprot_id = response["results"][0]["primaryAccession"]

            if uniprot_id:
                # 通过PDBe API获取结构映射
                url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
                response = requests.get(url)
                if response.status_code == 200:
                    structures = response.json().get(uniprot_id, [])
                    sorted_structures = sorted(structures, key=lambda x: x.get('resolution') or 999)
                    return [s['pdb_id'] for s in sorted_structures[:max_structures]]
        except Exception as e:
            print(f"⚠️  映射过程中出现错误: {e}")

        return []

    # ==================== 2. PDB查询与获取 ====================
    def fetch_pdb_info(self, pdb_id):
        """获取PDB结构详细信息"""
        url = f"{self.rcsb_base}/core/entry/{pdb_id}"
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                # 提取关键信息
                info = {
                    'pdb_id': pdb_id,
                    'title': data.get('struct', {}).get('title', 'N/A'),
                    'resolution': data.get('rcsb_entry_info', {}).get('resolution_combined', ['N/A'])[0],
                    'method': data.get('exptl', [{}])[0].get('method', 'N/A'),
                    'organism': data.get('rcsb_entity_source_organism', [{}])[0].get('scientific_name', 'N/A'),
                    'release_date': data.get('rcsb_accession_info', {}).get('deposit_date', 'N/A'),
                    'chains': []
                }

                # 获取链信息
                polymer_url = f"{self.rcsb_base}/core/polymer_entity/{pdb_id}/1"
                polymer_resp = requests.get(polymer_url)
                if polymer_resp.status_code == 200:
                    polymer_data = polymer_resp.json()
                    info['sequence'] = polymer_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
                    info['length'] = len(info['sequence']) if info['sequence'] else 0

                return info
        except Exception as e:
            print(f"获取PDB信息失败: {e}")
        return None

    # ==================== 3. 3D可视化与对比 ====================
    def view_3d(self, pdb_id, style='cartoon', color='spectrum', surface=False):
        """3D可视化单个结构"""
        viewer = py3Dmol.view(query=f'pdb:{pdb_id}')

        styles = {
            'cartoon': {'cartoon': {'color': color}},
            'stick': {'stick': {'colorscheme': 'greenCarbon'}},
            'sphere': {'sphere': {'radius': 0.5}}
        }

        viewer.setStyle(styles.get(style, styles['cartoon']))
        if surface:
            viewer.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'white'})

        viewer.zoomTo()
        print(f"✅ 正在加载 {pdb_id} 的3D结构...")
        return viewer

    def compare_structures(self, pdb_id1, pdb_id2, align=False):
        """对比两个结构"""
        viewer = py3Dmol.view()

        # 获取结构数据
        pdb_data1 = requests.get(f'https://files.rcsb.org/view/{pdb_id1}.pdb').text
        pdb_data2 = requests.get(f'https://files.rcsb.org/view/{pdb_id2}.pdb').text

        viewer.addModel(pdb_data1, 'pdb')
        viewer.setStyle({'model': 0}, {'cartoon': {'color': 'red'}})

        viewer.addModel(pdb_data2, 'pdb')
        viewer.setStyle({'model': 1}, {'cartoon': {'color': 'blue'}})

        if align:
            # 简单对齐（专业对齐需使用BioPython的Superimposer）
            viewer.align({'model': 1}, {'model': 0})

        viewer.zoomTo()
        print(f"🔄 正在对比 {pdb_id1} (红色) 和 {pdb_id2} (蓝色)")
        return viewer

    # ==================== 4. 物化性质分析 ====================
    def analyze_structure(self, pdb_id, properties=['all']):
        """分析蛋白结构的物化性质"""
        print(f"🧪 正在分析 {pdb_id} 的物化性质...")

        # 下载PDB文件
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

        if not pdb_file:
            return None

        # 解析结构
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]

        results = {'pdb_id': pdb_id}

        # 1. 基础信息
        results['num_chains'] = len(list(model.get_chains()))
        results['num_residues'] = len(list(model.get_residues()))
        results['num_atoms'] = len(list(model.get_atoms()))

        # 2. 序列分析（如果可用）
        if hasattr(self, 'sequence') and self.sequence:
            protein_analyzer = ProtParam.ProteinAnalysis(self.sequence)
            results['molecular_weight'] = protein_analyzer.molecular_weight()
            results['isoelectric_point'] = protein_analyzer.isoelectric_point()
            results['amino_acid_composition'] = protein_analyzer.get_amino_acids_percent()

        # 3. 二级结构估算（需要DSSP，此处为简化版）
        try:
            dssp = DSSP(model, pdb_file)
            ss_counts = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0}
            for key in dssp.keys():
                ss = dssp[key][2]
                if ss in ss_counts:
                    ss_counts[ss] += 1

            results['secondary_structure'] = {
                'helix': ss_counts['H'],
                'beta_sheet': ss_counts['E'],
                'coil': sum(ss_counts.values()) - (ss_counts['H'] + ss_counts['E'])
            }
        except:
            results['secondary_structure'] = {'helix': 'N/A', 'beta_sheet': 'N/A', 'coil': 'N/A'}

        return results

    # ==================== 5. 报告生成 ====================
    def generate_report(self, gene_name=None, pdb_ids=None):
        """生成交互式分析报告"""
        report = []
        report.append("# 🧬 蛋白结构综合分析报告\n")

        if gene_name:
            report.append(f"## 1. 基因查询: {gene_name}")
            structures = self.gene_to_structures(gene_name)

            if structures:
                report.append(f"找到 {len(structures)} 个相关结构:")
                for i, pdb_id in enumerate(structures[:3], 1):
                    info = self.fetch_pdb_info(pdb_id)
                    if info:
                        report.append(f"{i}. **{pdb_id}**: {info['title']} (分辨率: {info['resolution']}Å)")
                pdb_ids = structures[:2]  # 取前两个进行分析
            else:
                report.append("⚠️ 未找到相关结构，请直接提供PDB ID")
                pdb_ids = pdb_ids or []

        if pdb_ids:
            report.append("\n## 2. 结构分析")

            # 分析每个结构
            for i, pdb_id in enumerate(pdb_ids[:2]):  # 限制数量
                report.append(f"\n### 结构 {i + 1}: {pdb_id}")

                # 基本信息
                info = self.fetch_pdb_info(pdb_id)
                if info:
                    report.append(f"- **标题**: {info['title']}")
                    report.append(f"- **分辨率**: {info['resolution']}Å")
                    report.append(f"- **实验方法**: {info['method']}")
                    report.append(f"- **来源生物**: {info['organism']}")

                # 物化性质
                analysis = self.analyze_structure(pdb_id)
                if analysis:
                    report.append("\n**物化性质**:")
                    report.append(f"- 链数: {analysis['num_chains']}")
                    report.append(f"- 残基数: {analysis['num_residues']}")
                    report.append(f"- 原子数: {analysis['num_atoms']}")

        # 可视化部分
        report.append("\n## 3. 3D可视化")
        report.append("运行以下代码查看3D结构:")
        report.append("```python")
        if pdb_ids:
            report.append(f"# 查看单个结构\nviewer = gget_pdb.view_3d('{pdb_ids[0]}')")
            if len(pdb_ids) > 1:
                report.append(
                    f"\n# 对比两个结构\ncomparison = gget_pdb.compare_structures('{pdb_ids[0]}', '{pdb_ids[1]}')")
        report.append("```")

        # 在线链接
        report.append("\n## 4. 在线查看")
        if pdb_ids:
            for pdb_id in pdb_ids[:3]:
                report.append(f"- [{pdb_id} RCSB官方查看器](https://www.rcsb.org/3d-view/{pdb_id})")
                report.append(f"- [{pdb_id} Molstar查看器](https://molstar.org/viewer/?pdb-id={pdb_id})")

        return "\n".join(report)

    # ==================== 便捷函数 ====================
    def quick_analysis(self, input_term):
        """一键式快速分析：接受基因名或PDB ID"""
        result = {}

        # 判断输入类型
        if len(input_term) == 4 and input_term.isalnum():  # 可能是PDB ID
            result['type'] = 'pdb_id'
            result['pdb_ids'] = [input_term]
            result['info'] = self.fetch_pdb_info(input_term)
            result['analysis'] = self.analyze_structure(input_term)
        else:  # 可能是基因名
            result['type'] = 'gene'
            result['gene_name'] = input_term
            result['pdb_ids'] = self.gene_to_structures(input_term)
            if result['pdb_ids']:
                result['info'] = self.fetch_pdb_info(result['pdb_ids'][0])
                result['analysis'] = self.analyze_structure(result['pdb_ids'][0])

        return result


# 创建全局实例
gget_pdb = GGETPDB()


# 便捷函数别名
def pdb_view(pdb_id, **kwargs):
    return gget_pdb.view_3d(pdb_id, **kwargs)


def gene_view(gene_name, **kwargs):
    return gget_pdb.quick_analysis(gene_name)