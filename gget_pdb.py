#文件:gget_pdb.py
import requests
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.DSSP import DSSP
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import ProtParam
from Bio.Align import PairwiseAligner
import py3Dmol
import pandas as pd
import warnings
import re

warnings.filterwarnings('ignore')

# 氨基酸属性常量
AMINO_ACID_PROPERTIES = {
    'A': {'name': 'Alanine', 'charge': 0, 'hydrophobic': True, 'volume': 88.6, 'polar': False},
    'C': {'name': 'Cysteine', 'charge': 0, 'hydrophobic': True, 'volume': 108.5, 'polar': False},
    'D': {'name': 'Aspartic acid', 'charge': -1, 'hydrophobic': False, 'volume': 111.1, 'polar': True},
    'E': {'name': 'Glutamic acid', 'charge': -1, 'hydrophobic': False, 'volume': 138.4, 'polar': True},
    'F': {'name': 'Phenylalanine', 'charge': 0, 'hydrophobic': True, 'volume': 189.9, 'polar': False},
    'G': {'name': 'Glycine', 'charge': 0, 'hydrophobic': True, 'volume': 60.1, 'polar': False},
    'H': {'name': 'Histidine', 'charge': 0.5, 'hydrophobic': False, 'volume': 153.2, 'polar': True},
    'I': {'name': 'Isoleucine', 'charge': 0, 'hydrophobic': True, 'volume': 166.7, 'polar': False},
    'K': {'name': 'Lysine', 'charge': 1, 'hydrophobic': False, 'volume': 168.6, 'polar': True},
    'L': {'name': 'Leucine', 'charge': 0, 'hydrophobic': True, 'volume': 166.7, 'polar': False},
    'M': {'name': 'Methionine', 'charge': 0, 'hydrophobic': True, 'volume': 162.9, 'polar': False},
    'N': {'name': 'Asparagine', 'charge': 0, 'hydrophobic': False, 'volume': 114.1, 'polar': True},
    'P': {'name': 'Proline', 'charge': 0, 'hydrophobic': True, 'volume': 112.7, 'polar': False},
    'Q': {'name': 'Glutamine', 'charge': 0, 'hydrophobic': False, 'volume': 143.8, 'polar': True},
    'R': {'name': 'Arginine', 'charge': 1, 'hydrophobic': False, 'volume': 173.4, 'polar': True},
    'S': {'name': 'Serine', 'charge': 0, 'hydrophobic': False, 'volume': 89.0, 'polar': True},
    'T': {'name': 'Threonine', 'charge': 0, 'hydrophobic': False, 'volume': 116.1, 'polar': True},
    'V': {'name': 'Valine', 'charge': 0, 'hydrophobic': True, 'volume': 140.0, 'polar': False},
    'W': {'name': 'Tryptophan', 'charge': 0, 'hydrophobic': True, 'volume': 227.8, 'polar': False},
    'Y': {'name': 'Tyrosine', 'charge': 0, 'hydrophobic': False, 'volume': 193.6, 'polar': True},
}

# 三字母到单字母氨基酸转换
THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}


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
                    'organism': 'N/A',  # 将从polymer_entity获取
                    'release_date': data.get('rcsb_accession_info', {}).get('deposit_date', 'N/A'),
                    'chains': []
                }

                # 获取链信息和来源生物
                polymer_url = f"{self.rcsb_base}/core/polymer_entity/{pdb_id}/1"
                polymer_resp = requests.get(polymer_url)
                if polymer_resp.status_code == 200:
                    polymer_data = polymer_resp.json()
                    info['sequence'] = polymer_data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
                    info['length'] = len(info['sequence']) if info['sequence'] else 0

                    # 从polymer_entity获取来源生物信息
                    # 优先使用 rcsb_entity_source_organism，如果没有则使用 rcsb_entity_host_organism
                    source_organism = polymer_data.get('rcsb_entity_source_organism')
                    if source_organism and len(source_organism) > 0:
                        info['organism'] = source_organism[0].get('scientific_name', 'N/A')
                    else:
                        host_organism = polymer_data.get('rcsb_entity_host_organism')
                        if host_organism and len(host_organism) > 0:
                            info['organism'] = host_organism[0].get('scientific_name', 'N/A')

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
    def analyze_structure(self, pdb_id, properties=None):
        """分析蛋白结构的物化性质"""
        if properties is None:
            properties = ['all']
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

        # 3. 二级结构估算（优先从API获取，失败后尝试DSSP）
        # 首先尝试从API获取二级结构
        ss_from_api = self._get_secondary_structure_from_api(pdb_id)
        if ss_from_api:
            results['secondary_structure'] = ss_from_api
        else:
            # 备用方案：尝试DSSP
            try:
                dssp = DSSP(model, pdb_file)
                ss_counts = {'H': 0, 'B': 0, 'E': 0, 'G': 0, 'I': 0, 'T': 0, 'S': 0, '-': 0}
                for key in dssp.keys():
                    ss = dssp[key][2]
                    if ss in ss_counts:
                        ss_counts[ss] += 1

                total = sum(ss_counts.values())
                results['secondary_structure'] = {
                    'helix': ss_counts['H'] + ss_counts['G'] + ss_counts['I'],
                    'beta_sheet': ss_counts['E'] + ss_counts['B'],
                    'coil': ss_counts['T'] + ss_counts['S'] + ss_counts['-'],
                    'helix_pct': round((ss_counts['H'] + ss_counts['G'] + ss_counts['I']) / total * 100, 1) if total > 0 else 0,
                    'beta_pct': round((ss_counts['E'] + ss_counts['B']) / total * 100, 1) if total > 0 else 0,
                    'coil_pct': round((ss_counts['T'] + ss_counts['S'] + ss_counts['-']) / total * 100, 1) if total > 0 else 0,
                    'source': 'DSSP'
                }
            except Exception as e:
                results['secondary_structure'] = {
                    'helix': 'N/A',
                    'beta_sheet': 'N/A',
                    'coil': 'N/A',
                    'note': 'DSSP未安装，请运行: brew install dssp (macOS) 或 apt-get install dssp (Linux)'
                }

        return results

    def _get_secondary_structure_from_api(self, pdb_id):
        """从RCSB/PDBe API获取二级结构信息"""
        # 方法1: 尝试从PDBe API获取二级结构注解
        try:
            pdbe_url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure/{pdb_id.lower()}"
            response = requests.get(pdbe_url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                pdb_data = data.get(pdb_id.lower(), {})
                molecules = pdb_data.get('molecules', [])

                total_helix = 0
                total_strand = 0
                total_coil = 0
                total_residues = 0

                for molecule in molecules:
                    chains = molecule.get('chains', [])
                    for chain in chains:
                        secondary_structure = chain.get('secondary_structure', {})

                        # 统计螺旋
                        helices = secondary_structure.get('helices', [])
                        for helix in helices:
                            start = helix.get('start', {}).get('residue_number', 0)
                            end = helix.get('end', {}).get('residue_number', 0)
                            if end >= start:
                                total_helix += (end - start + 1)

                        # 统计β折叠
                        strands = secondary_structure.get('strands', [])
                        for strand in strands:
                            start = strand.get('start', {}).get('residue_number', 0)
                            end = strand.get('end', {}).get('residue_number', 0)
                            if end >= start:
                                total_strand += (end - start + 1)

                # 从RCSB获取总残基数
                try:
                    rcsb_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                    rcsb_resp = requests.get(rcsb_url, timeout=10)
                    if rcsb_resp.status_code == 200:
                        rcsb_data = rcsb_resp.json()
                        total_residues = rcsb_data.get('rcsb_entry_info', {}).get('deposited_polymer_monomer_count', 0)
                except:
                    pass

                if total_helix > 0 or total_strand > 0:
                    total_coil = max(0, total_residues - total_helix - total_strand) if total_residues > 0 else 0

                    result = {
                        'helix': total_helix,
                        'beta_sheet': total_strand,
                        'coil': total_coil if total_residues > 0 else 'N/A',
                        'source': 'PDBe API'
                    }

                    # 计算百分比
                    if total_residues > 0:
                        result['helix_pct'] = round(total_helix / total_residues * 100, 1)
                        result['beta_pct'] = round(total_strand / total_residues * 100, 1)
                        result['coil_pct'] = round(total_coil / total_residues * 100, 1)

                    return result
        except Exception as e:
            print(f"PDBe API获取二级结构失败: {e}")

        # 方法2: 尝试从RCSB获取简化的二级结构信息
        try:
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/1"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                # 从entity_poly获取序列长度
                seq_length = len(data.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', ''))

                # 尝试获取二级结构注解
                annotations = data.get('rcsb_polymer_entity_annotation', [])
                helix_residues = 0
                sheet_residues = 0

                for annotation in annotations:
                    ann_type = annotation.get('type', '')
                    if 'HELIX' in ann_type.upper():
                        # 尝试获取范围
                        feature = annotation.get('annotation_lineage', [{}])
                        helix_residues += 1
                    elif 'SHEET' in ann_type.upper() or 'STRAND' in ann_type.upper():
                        sheet_residues += 1

                if helix_residues > 0 or sheet_residues > 0:
                    return {
                        'helix': helix_residues,
                        'beta_sheet': sheet_residues,
                        'coil': 'N/A',
                        'source': 'RCSB API (简化)'
                    }
        except Exception as e:
            print(f"RCSB API获取二级结构失败: {e}")

        return None

    def _get_hydrogen_bonds_from_api(self, pdb_id):
        """从PDBe API获取氢键信息"""
        try:
            # PDBe提供的分子间相互作用API
            url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id.lower()}"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                pdb_data = data.get(pdb_id.lower(), [{}])[0]

                # 从摘要中尝试获取相关信息
                # 注意：PDBe API不直接提供氢键数量，这里返回估算信息
                num_entities = pdb_data.get('number_of_entities', {}).get('polypeptide', 0)
                total_residues = 0

                # 获取残基数用于估算
                try:
                    rcsb_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
                    rcsb_resp = requests.get(rcsb_url, timeout=10)
                    if rcsb_resp.status_code == 200:
                        rcsb_data = rcsb_resp.json()
                        total_residues = rcsb_data.get('rcsb_entry_info', {}).get('deposited_polymer_monomer_count', 0)
                except:
                    pass

                if total_residues > 0:
                    # 粗略估算：平均每个残基约有0.7-1.0个主链氢键
                    estimated_hbonds = int(total_residues * 0.85)
                    return {
                        'backbone_hbonds': f"~{estimated_hbonds} (估算)",
                        'total': f"~{estimated_hbonds} (估算)",
                        'source': 'PDBe API (估算值)',
                        'note': '基于残基数估算，精确值需要DSSP分析'
                    }
        except Exception as e:
            print(f"PDBe API获取氢键信息失败: {e}")

        return None

    # ==================== 4.1 高级结构分析 ====================
    def analyze_advanced_structure(self, pdb_id):
        """高级结构分析：氢键、盐桥、二硫键、SASA、疏水/亲水比例"""
        print(f"🔬 正在进行 {pdb_id} 的高级结构分析...")

        # 下载PDB文件
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

        if not pdb_file:
            return None

        # 解析结构
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]

        results = {'pdb_id': pdb_id}

        # 1. 二硫键分析
        results['disulfide_bonds'] = self._find_disulfide_bonds(model)

        # 2. 盐桥分析
        results['salt_bridges'] = self._find_salt_bridges(model)

        # 3. 氢键统计（优先从API获取，失败后尝试DSSP）
        # 首先尝试从API获取氢键信息
        hbonds_from_api = self._get_hydrogen_bonds_from_api(pdb_id)
        if hbonds_from_api:
            results['hydrogen_bonds'] = hbonds_from_api
        else:
            # 备用方案：尝试DSSP
            try:
                dssp = DSSP(model, pdb_file)
                results['hydrogen_bonds'] = self._count_hydrogen_bonds(dssp)
                results['hydrogen_bonds']['source'] = 'DSSP'
            except Exception as e:
                # 提供更详细的错误信息和解决方案
                error_msg = str(e)
                if 'mkdssp' in error_msg.lower() or 'dssp' in error_msg.lower():
                    results['hydrogen_bonds'] = {
                        'backbone_hbonds': 'N/A',
                        'total': 'N/A',
                        'note': 'DSSP未安装，请运行: brew install dssp (macOS) 或 apt-get install dssp (Linux)'
                    }
                else:
                    results['hydrogen_bonds'] = {
                        'backbone_hbonds': 'N/A',
                        'total': 'N/A',
                        'error': error_msg
                    }

        # 4. SASA分析（每条链）
        results['sasa_per_chain'] = self._calculate_sasa(model)

        # 5. 疏水/亲水残基比例（每条链）
        results['hydrophobicity_per_chain'] = self._analyze_hydrophobicity(model)

        return results

    def _find_disulfide_bonds(self, model):
        """查找二硫键"""
        disulfide_bonds = []
        cysteine_residues = []

        # 收集所有半胱氨酸的SG原子
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    if 'SG' in residue:
                        cysteine_residues.append({
                            'chain': chain.id,
                            'resnum': residue.id[1],
                            'atom': residue['SG']
                        })

        # 检查半胱氨酸之间的距离（二硫键距离约2.05Å）
        for i, cys1 in enumerate(cysteine_residues):
            for cys2 in cysteine_residues[i+1:]:
                distance = float(cys1['atom'] - cys2['atom'])  # 转换为 Python float
                if distance < 2.5:  # 二硫键距离阈值
                    disulfide_bonds.append({
                        'cys1': f"{cys1['chain']}:{cys1['resnum']}",
                        'cys2': f"{cys2['chain']}:{cys2['resnum']}",
                        'distance': round(distance, 2)
                    })

        return {'count': len(disulfide_bonds), 'bonds': disulfide_bonds}

    def _find_salt_bridges(self, model, distance_cutoff=4.0):
        """查找盐桥"""
        salt_bridges = []

        # 正电荷残基的原子
        positive_atoms = []
        # 负电荷残基的原子
        negative_atoms = []

        positive_residues = ['ARG', 'LYS', 'HIS']
        negative_residues = ['ASP', 'GLU']

        positive_atom_names = {'ARG': ['NH1', 'NH2', 'NE'], 'LYS': ['NZ'], 'HIS': ['ND1', 'NE2']}
        negative_atom_names = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}

        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if resname in positive_residues:
                    for atom_name in positive_atom_names.get(resname, []):
                        if atom_name in residue:
                            positive_atoms.append({
                                'chain': chain.id,
                                'resname': resname,
                                'resnum': residue.id[1],
                                'atom': residue[atom_name]
                            })
                elif resname in negative_residues:
                    for atom_name in negative_atom_names.get(resname, []):
                        if atom_name in residue:
                            negative_atoms.append({
                                'chain': chain.id,
                                'resname': resname,
                                'resnum': residue.id[1],
                                'atom': residue[atom_name]
                            })

        # 检查正负电荷原子之间的距离
        seen_pairs = set()
        for pos in positive_atoms:
            for neg in negative_atoms:
                distance = float(pos['atom'] - neg['atom'])  # 转换为 Python float
                if distance <= distance_cutoff:
                    pair_key = tuple(sorted([
                        f"{pos['chain']}:{pos['resname']}{pos['resnum']}",
                        f"{neg['chain']}:{neg['resname']}{neg['resnum']}"
                    ]))
                    if pair_key not in seen_pairs:
                        seen_pairs.add(pair_key)
                        salt_bridges.append({
                            'positive': f"{pos['chain']}:{pos['resname']}{pos['resnum']}",
                            'negative': f"{neg['chain']}:{neg['resname']}{neg['resnum']}",
                            'distance': round(distance, 2)
                        })

        return {'count': len(salt_bridges), 'bridges': salt_bridges}

    def _count_hydrogen_bonds(self, dssp):
        """统计氢键数量（基于DSSP）"""
        # DSSP提供的氢键信息
        hbonds_backbone = 0

        for key in dssp.keys():
            # DSSP返回的氢键信息（NH-->O和O-->NH）
            # 索引3和4是NH-->O氢键，5和6是O-->NH氢键
            dssp_data = dssp[key]
            # 检查NH-->O方向
            if dssp_data[6] != 0:  # 能量不为0表示存在氢键
                hbonds_backbone += 1
            if dssp_data[8] != 0:
                hbonds_backbone += 1

        return {'backbone_hbonds': hbonds_backbone, 'total': hbonds_backbone}

    def _calculate_sasa(self, model):
        """计算每条链的SASA"""
        sasa_results = {}

        try:
            # 使用ShrakeRupley算法计算SASA
            sr = ShrakeRupley()
            sr.compute(model, level="R")  # 残基级别

            for chain in model:
                chain_id = chain.id
                total_sasa = 0.0
                for residue in chain:
                    if hasattr(residue, 'sasa'):
                        total_sasa += float(residue.sasa)  # 转换为 Python float

                sasa_results[chain_id] = round(float(total_sasa), 2)  # 确保是 Python float
        except Exception as e:
            sasa_results['error'] = str(e)

        return sasa_results

    def _analyze_hydrophobicity(self, model):
        """分析每条链的疏水/亲水残基比例"""
        results = {}

        for chain in model:
            chain_id = chain.id
            hydrophobic_count = 0
            hydrophilic_count = 0
            total = 0

            for residue in chain:
                resname = residue.get_resname()
                one_letter = THREE_TO_ONE.get(resname)
                if one_letter and one_letter in AMINO_ACID_PROPERTIES:
                    total += 1
                    if AMINO_ACID_PROPERTIES[one_letter]['hydrophobic']:
                        hydrophobic_count += 1
                    else:
                        hydrophilic_count += 1

            if total > 0:
                results[chain_id] = {
                    'hydrophobic_count': hydrophobic_count,
                    'hydrophilic_count': hydrophilic_count,
                    'hydrophobic_ratio': round(hydrophobic_count / total * 100, 2),
                    'hydrophilic_ratio': round(hydrophilic_count / total * 100, 2),
                    'total_residues': total
                }

        return results

    # ==================== 4.2 突变影响分析 ====================
    def analyze_mutation(self, pdb_id, mutation_str):
        """
        分析突变影响
        mutation_str格式: "A:K33E" 表示A链第33位由K突变为E
        """
        print(f"🧬 正在分析突变 {mutation_str} 对 {pdb_id} 的影响...")

        # 解析突变字符串
        match = re.match(r'([A-Z]):([A-Z])(\d+)([A-Z])', mutation_str.upper())
        if not match:
            return {'error': '突变格式无效，请使用格式: A:K33E (链:原氨基酸+位置+新氨基酸)'}

        chain_id, wt_aa, position, mut_aa = match.groups()
        position = int(position)

        # 验证氨基酸
        if wt_aa not in AMINO_ACID_PROPERTIES or mut_aa not in AMINO_ACID_PROPERTIES:
            return {'error': f'无效的氨基酸代码: {wt_aa} 或 {mut_aa}'}

        wt_props = AMINO_ACID_PROPERTIES[wt_aa]
        mut_props = AMINO_ACID_PROPERTIES[mut_aa]

        # 计算变化
        charge_change = mut_props['charge'] - wt_props['charge']
        volume_change = mut_props['volume'] - wt_props['volume']
        hydrophobicity_change = mut_props['hydrophobic'] != wt_props['hydrophobic']
        polarity_change = mut_props['polar'] != wt_props['polar']

        # 评估影响
        impact_score = 0
        impact_reasons = []

        if abs(charge_change) >= 1:
            impact_score += 3
            impact_reasons.append(f"电荷变化: {'+' if charge_change > 0 else ''}{charge_change}")

        if abs(volume_change) > 50:
            impact_score += 2
            impact_reasons.append(f"体积变化: {'+' if volume_change > 0 else ''}{volume_change:.1f}Å³")
        elif abs(volume_change) > 20:
            impact_score += 1
            impact_reasons.append(f"中等体积变化: {'+' if volume_change > 0 else ''}{volume_change:.1f}Å³")

        if hydrophobicity_change:
            impact_score += 2
            if wt_props['hydrophobic']:
                impact_reasons.append("疏水 → 亲水 (可能影响蛋白折叠)")
            else:
                impact_reasons.append("亲水 → 疏水 (可能影响溶解性)")

        if polarity_change:
            impact_score += 1
            impact_reasons.append("极性变化")

        # 下载并检查结构中的实际残基
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

        structural_context = None
        if pdb_file:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(pdb_id, pdb_file)
            model = structure[0]

            try:
                chain = model[chain_id]
                residue = chain[position]
                actual_resname = THREE_TO_ONE.get(residue.get_resname(), '?')

                structural_context = {
                    'found_residue': actual_resname,
                    'matches_wt': actual_resname == wt_aa,
                    'position_valid': True
                }

                if actual_resname != wt_aa:
                    structural_context['warning'] = f"结构中该位置的氨基酸是 {actual_resname}，而非 {wt_aa}"

                # 检查是否在二级结构中
                try:
                    dssp = DSSP(model, pdb_file)
                    dssp_key = (chain_id, (' ', position, ' '))
                    if dssp_key in dssp:
                        ss = dssp[dssp_key][2]
                        ss_mapping = {
                            'H': 'α-螺旋', 'G': '3₁₀-螺旋', 'I': 'π-螺旋',
                            'E': 'β-折叠', 'B': 'β-桥', 'T': '转角',
                            'S': '弯曲', '-': '环区'
                        }
                        structural_context['secondary_structure'] = ss_mapping.get(ss, ss)

                        # 在二级结构核心区域的突变影响更大
                        if ss in ['H', 'E']:
                            impact_score += 1
                            impact_reasons.append(f"位于{ss_mapping[ss]}核心区域")
                except:
                    pass

            except KeyError:
                structural_context = {
                    'found_residue': None,
                    'matches_wt': False,
                    'position_valid': False,
                    'error': f"未找到链 {chain_id} 或位置 {position}"
                }

        # 生成影响评估
        if impact_score >= 5:
            impact_level = "高"
            impact_description = "该突变可能严重影响蛋白结构或功能"
        elif impact_score >= 3:
            impact_level = "中"
            impact_description = "该突变可能对蛋白有中等程度的影响"
        else:
            impact_level = "低"
            impact_description = "该突变可能是保守性替换，影响较小"

        return {
            'mutation': mutation_str,
            'pdb_id': pdb_id,
            'wild_type': {
                'aa': wt_aa,
                'name': wt_props['name'],
                'charge': wt_props['charge'],
                'volume': wt_props['volume'],
                'hydrophobic': wt_props['hydrophobic']
            },
            'mutant': {
                'aa': mut_aa,
                'name': mut_props['name'],
                'charge': mut_props['charge'],
                'volume': mut_props['volume'],
                'hydrophobic': mut_props['hydrophobic']
            },
            'changes': {
                'charge_change': charge_change,
                'volume_change': round(volume_change, 2),
                'hydrophobicity_change': hydrophobicity_change,
                'polarity_change': polarity_change
            },
            'impact_assessment': {
                'score': impact_score,
                'level': impact_level,
                'description': impact_description,
                'reasons': impact_reasons
            },
            'structural_context': structural_context
        }

    # ==================== 4.3 序列分析 ====================
    def analyze_sequence_composition(self, pdb_id):
        """分析每条链的氨基酸组成"""
        print(f"📊 正在分析 {pdb_id} 的序列组成...")

        # 下载PDB文件
        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

        if not pdb_file:
            return None

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]

        results = {'pdb_id': pdb_id, 'chains': {}}

        for chain in model:
            chain_id = chain.id
            sequence = []
            aa_counts = {aa: 0 for aa in AMINO_ACID_PROPERTIES.keys()}

            for residue in chain:
                resname = residue.get_resname()
                one_letter = THREE_TO_ONE.get(resname)
                if one_letter:
                    sequence.append(one_letter)
                    if one_letter in aa_counts:
                        aa_counts[one_letter] += 1

            if sequence:
                total = len(sequence)
                # 计算百分比
                aa_percentages = {aa: round(count / total * 100, 2)
                                 for aa, count in aa_counts.items()}

                # 分类统计
                charged_positive = sum(aa_counts[aa] for aa in ['K', 'R', 'H'])
                charged_negative = sum(aa_counts[aa] for aa in ['D', 'E'])
                hydrophobic = sum(aa_counts[aa] for aa in ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'])
                polar = sum(aa_counts[aa] for aa in ['S', 'T', 'N', 'Q', 'Y', 'C'])
                aromatic = sum(aa_counts[aa] for aa in ['F', 'Y', 'W'])

                results['chains'][chain_id] = {
                    'sequence': ''.join(sequence),
                    'length': total,
                    'amino_acid_counts': aa_counts,
                    'amino_acid_percentages': aa_percentages,
                    'category_statistics': {
                        'charged_positive': charged_positive,
                        'charged_positive_pct': round(charged_positive / total * 100, 2),
                        'charged_negative': charged_negative,
                        'charged_negative_pct': round(charged_negative / total * 100, 2),
                        'hydrophobic': hydrophobic,
                        'hydrophobic_pct': round(hydrophobic / total * 100, 2),
                        'polar_uncharged': polar,
                        'polar_uncharged_pct': round(polar / total * 100, 2),
                        'aromatic': aromatic,
                        'aromatic_pct': round(aromatic / total * 100, 2)
                    }
                }

        return results

    def align_with_uniprot(self, pdb_id, uniprot_id=None):
        """将PDB序列与UniProt canonical序列比对"""
        print(f"🔗 正在比对 {pdb_id} 与 UniProt 序列...")

        # 获取PDB序列
        pdb_info = self.fetch_pdb_info(pdb_id)
        if not pdb_info or not pdb_info.get('sequence'):
            # 尝试从结构文件获取
            pdbl = PDBList()
            pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
            if pdb_file:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure(pdb_id, pdb_file)
                model = structure[0]

                pdb_sequences = {}
                for chain in model:
                    seq = []
                    for residue in chain:
                        one_letter = THREE_TO_ONE.get(residue.get_resname())
                        if one_letter:
                            seq.append(one_letter)
                    if seq:
                        pdb_sequences[chain.id] = ''.join(seq)
            else:
                return {'error': f'无法获取 {pdb_id} 的序列'}
        else:
            pdb_sequences = {'A': pdb_info.get('sequence', '')}

        # 如果没有提供UniProt ID，尝试从PDB映射获取
        if not uniprot_id:
            try:
                url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
                response = requests.get(url)
                if response.status_code == 200:
                    data = response.json()
                    uniprot_entries = data.get(pdb_id.lower(), {}).get('UniProt', {})
                    if uniprot_entries:
                        uniprot_id = list(uniprot_entries.keys())[0]
            except:
                pass

        if not uniprot_id:
            return {'error': '无法确定UniProt ID，请手动提供', 'pdb_sequences': pdb_sequences}

        # 获取UniProt序列
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            response = requests.get(url)
            if response.status_code != 200:
                return {'error': f'无法获取UniProt序列: {uniprot_id}'}

            fasta_lines = response.text.strip().split('\n')
            uniprot_seq = ''.join(fasta_lines[1:])
        except Exception as e:
            return {'error': f'获取UniProt序列失败: {e}'}

        # 进行序列比对
        results = {
            'pdb_id': pdb_id,
            'uniprot_id': uniprot_id,
            'uniprot_length': len(uniprot_seq),
            'chain_alignments': {}
        }

        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        for chain_id, pdb_seq in pdb_sequences.items():
            if len(pdb_seq) < 10:  # 跳过太短的序列
                continue

            alignments = aligner.align(uniprot_seq, pdb_seq)
            if alignments:
                best_alignment = alignments[0]

                # 计算序列一致性
                aligned_uniprot, aligned_pdb = best_alignment.aligned

                # 简单计算一致性
                matches = 0
                total_aligned = 0
                gaps_in_pdb = []  # 缺失区段
                insertions_in_pdb = []  # 插入区段

                uniprot_aligned = str(best_alignment).split('\n')[0]
                pdb_aligned = str(best_alignment).split('\n')[2] if len(str(best_alignment).split('\n')) > 2 else ''

                # 计算identity
                for i, (u, p) in enumerate(zip(uniprot_aligned, pdb_aligned)):
                    if u != '-' and p != '-':
                        total_aligned += 1
                        if u == p:
                            matches += 1

                identity = round(matches / len(uniprot_seq) * 100, 2) if uniprot_seq else 0
                coverage = round(len(pdb_seq) / len(uniprot_seq) * 100, 2) if uniprot_seq else 0

                # 检测缺失和插入区段
                # 通过aligned块来识别
                current_pos = 0
                for block in aligned_uniprot:
                    start, end = block
                    if start > current_pos:
                        gaps_in_pdb.append({'start': current_pos + 1, 'end': start, 'length': start - current_pos})
                    current_pos = end

                if current_pos < len(uniprot_seq):
                    gaps_in_pdb.append({'start': current_pos + 1, 'end': len(uniprot_seq), 'length': len(uniprot_seq) - current_pos})

                results['chain_alignments'][chain_id] = {
                    'pdb_length': len(pdb_seq),
                    'identity_percent': identity,
                    'coverage_percent': coverage,
                    'missing_regions': gaps_in_pdb,
                    'alignment_score': best_alignment.score
                }

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