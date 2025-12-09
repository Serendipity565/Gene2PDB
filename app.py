# 文件：app.py
# Flask 后端服务 API
from flask import Flask, jsonify, request
from flask_cors import CORS
from gget_pdb import GGETPDB

app = Flask(__name__)
CORS(app)  # 允许跨域请求

# 创建分析工具实例
analyzer = GGETPDB()


@app.route('/api/health', methods=['GET'])
def health_check():
    """健康检查接口"""
    return jsonify({'status': 'ok', 'message': 'PDB分析服务正常运行'})


@app.route('/api/gene/structures', methods=['GET'])
def get_gene_structures():
    """根据基因名查找相关PDB结构"""
    gene_name = request.args.get('gene_name', '')
    species = request.args.get('species', 'human')
    max_structures = int(request.args.get('max_structures', 5))

    if not gene_name:
        return jsonify({'error': '请提供基因名称'}), 400

    try:
        structures = analyzer.gene_to_structures(gene_name, species=species, max_structures=max_structures)
        return jsonify({
            'gene_name': gene_name,
            'species': species,
            'structures': structures,
            'count': len(structures)
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/pdb/info/<pdb_id>', methods=['GET'])
def get_pdb_info(pdb_id):
    """获取PDB结构详细信息"""
    try:
        info = analyzer.fetch_pdb_info(pdb_id)
        if info:
            return jsonify(info)
        else:
            return jsonify({'error': f'未找到PDB结构 {pdb_id}'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/pdb/analyze/<pdb_id>', methods=['GET'])
def analyze_pdb(pdb_id):
    """分析PDB结构的物化性质"""
    try:
        analysis = analyzer.analyze_structure(pdb_id)
        if analysis:
            return jsonify(analysis)
        else:
            return jsonify({'error': f'无法分析PDB结构 {pdb_id}'}), 404
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/report', methods=['GET'])
def generate_report():
    """生成分析报告"""
    gene_name = request.args.get('gene_name', '')
    pdb_ids = request.args.getlist('pdb_ids')

    try:
        report = analyzer.generate_report(gene_name=gene_name if gene_name else None,
                                         pdb_ids=pdb_ids if pdb_ids else None)
        return jsonify({'report': report})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/quick', methods=['GET'])
def quick_analysis():
    """快速分析：接受基因名或PDB ID"""
    input_term = request.args.get('input', '')

    if not input_term:
        return jsonify({'error': '请提供基因名或PDB ID'}), 400

    try:
        result = analyzer.quick_analysis(input_term)
        # 转换为可JSON序列化的格式
        serializable_result = {
            'type': result.get('type'),
            'gene_name': result.get('gene_name'),
            'pdb_ids': result.get('pdb_ids', []),
            'info': result.get('info'),
            'analysis': {
                k: v for k, v in (result.get('analysis') or {}).items()
                if not isinstance(v, (bytes, type(None))) or v is None
            } if result.get('analysis') else None
        }
        return jsonify(serializable_result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500


if __name__ == '__main__':
    print("🚀 PDB分析后端服务启动中...")
    print("📡 API文档:")
    print("   GET /api/health - 健康检查")
    print("   GET /api/gene/structures?gene_name=INS - 查找基因相关结构")
    print("   GET /api/pdb/info/<pdb_id> - 获取PDB信息")
    print("   GET /api/pdb/analyze/<pdb_id> - 分析PDB结构")
    print("   GET /api/report?gene_name=INS - 生成报告")
    print("   GET /api/quick?input=INS - 快速分析")
    app.run(debug=True, host='0.0.0.0', port=8080)

