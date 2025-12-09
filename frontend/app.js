// API 基础地址
const API_BASE = 'http://localhost:8080/api';

// 全局变量
let currentPdbId = null;
let viewer = null;
let currentReport = '';

// 页面加载完成后初始化
document.addEventListener('DOMContentLoaded', () => {
    // 绑定回车键搜索
    document.getElementById('searchInput').addEventListener('keypress', (e) => {
        if (e.key === 'Enter') {
            quickAnalysis();
        }
    });

    // 检查后端服务状态
    checkHealth();
});

// 检查后端服务健康状态
async function checkHealth() {
    try {
        const response = await fetch(`${API_BASE}/health`);
        const data = await response.json();
        console.log('后端服务状态:', data);
    } catch (error) {
        showError('无法连接到后端服务，请确保已启动 Flask 服务器 (python app.py)');
    }
}

// 快速分析
async function quickAnalysis() {
    const input = document.getElementById('searchInput').value.trim();
    if (!input) {
        showError('请输入基因名或 PDB ID');
        return;
    }

    showLoading(true);
    hideResults();

    const searchType = document.querySelector('input[name="searchType"]:checked').value;

    try {
        if (searchType === 'gene') {
            await searchByGene(input);
        } else {
            await searchByPdbId(input);
        }
    } catch (error) {
        showError(`分析失败: ${error.message}`);
    } finally {
        showLoading(false);
    }
}

// 按基因名搜索
async function searchByGene(geneName) {
    const species = document.getElementById('speciesSelect').value;

    // 获取相关结构
    const structuresResponse = await fetch(
        `${API_BASE}/gene/structures?gene_name=${encodeURIComponent(geneName)}&species=${species}`
    );
    const structuresData = await structuresResponse.json();

    if (structuresData.error) {
        throw new Error(structuresData.error);
    }

    if (!structuresData.structures || structuresData.structures.length === 0) {
        showError(`未找到基因 "${geneName}" 的相关 PDB 结构`);
        return;
    }

    // 显示基本信息
    displayBasicInfo({
        type: '基因搜索',
        gene_name: geneName,
        species: species,
        structure_count: structuresData.count
    });

    // 显示结构列表
    await displayStructuresList(structuresData.structures);

    // 生成报告
    await generateReport(geneName);

    showResults();
}

// 按 PDB ID 搜索
async function searchByPdbId(pdbId) {
    // 获取 PDB 信息
    const infoResponse = await fetch(`${API_BASE}/pdb/info/${pdbId}`);
    const infoData = await infoResponse.json();

    if (infoData.error) {
        throw new Error(infoData.error);
    }

    // 显示基本信息
    displayBasicInfo({
        type: 'PDB 搜索',
        pdb_id: pdbId,
        title: infoData.title,
        method: infoData.method
    });

    // 显示单个结构
    await displayStructuresList([pdbId]);

    // 显示详情
    await showStructureDetail(pdbId);

    // 生成报告
    await generateReportForPdb(pdbId);

    showResults();
}

// 显示基本信息
function displayBasicInfo(info) {
    const container = document.getElementById('basicInfo');
    let html = '';

    const labels = {
        type: '搜索类型',
        gene_name: '基因名',
        species: '物种',
        structure_count: '结构数量',
        pdb_id: 'PDB ID',
        title: '标题',
        method: '实验方法'
    };

    for (const [key, value] of Object.entries(info)) {
        if (value !== undefined && value !== null) {
            html += `
                <div class="info-item">
                    <div class="label">${labels[key] || key}</div>
                    <div class="value">${value}</div>
                </div>
            `;
        }
    }

    container.innerHTML = html;
}

// 显示结构列表
async function displayStructuresList(pdbIds) {
    const container = document.getElementById('structuresList');
    container.innerHTML = '<p>正在加载结构信息...</p>';

    let html = '';

    for (const pdbId of pdbIds) {
        try {
            const response = await fetch(`${API_BASE}/pdb/info/${pdbId}`);
            const info = await response.json();

            html += `
                <div class="structure-item" onclick="selectStructure('${pdbId}')" id="structure-${pdbId}">
                    <div class="pdb-id">${pdbId.toUpperCase()}</div>
                    <div class="title">${info.title || '无标题'}</div>
                    <div class="meta">
                        分辨率: ${info.resolution || 'N/A'}Å | 
                        方法: ${info.method || 'N/A'}
                    </div>
                </div>
            `;
        } catch (error) {
            html += `
                <div class="structure-item" onclick="selectStructure('${pdbId}')" id="structure-${pdbId}">
                    <div class="pdb-id">${pdbId.toUpperCase()}</div>
                    <div class="title">加载信息失败</div>
                </div>
            `;
        }
    }

    container.innerHTML = html;

    // 默认选择第一个结构
    if (pdbIds.length > 0) {
        selectStructure(pdbIds[0]);
    }
}

// 选择结构
async function selectStructure(pdbId) {
    // 更新选中状态
    document.querySelectorAll('.structure-item').forEach(item => {
        item.classList.remove('active');
    });
    const selectedItem = document.getElementById(`structure-${pdbId}`);
    if (selectedItem) {
        selectedItem.classList.add('active');
    }

    currentPdbId = pdbId;

    // 显示详情
    await showStructureDetail(pdbId);

    // 加载 3D 视图
    load3DViewer(pdbId);
}

// 显示结构详情
async function showStructureDetail(pdbId) {
    const container = document.getElementById('structureDetail');
    const content = document.getElementById('detailContent');

    container.classList.remove('hidden');
    content.innerHTML = '<p>正在加载详情...</p>';

    try {
        // 获取 PDB 信息
        const infoResponse = await fetch(`${API_BASE}/pdb/info/${pdbId}`);
        const info = await infoResponse.json();

        let html = `
            <div class="detail-section">
                <h3>📄 基本信息</h3>
                <p><strong>PDB ID:</strong> ${pdbId.toUpperCase()}</p>
                <p><strong>标题:</strong> ${info.title || 'N/A'}</p>
                <p><strong>分辨率:</strong> ${info.resolution || 'N/A'}Å</p>
                <p><strong>实验方法:</strong> ${info.method || 'N/A'}</p>
                <p><strong>来源生物:</strong> ${info.organism || 'N/A'}</p>
                <p><strong>发布日期:</strong> ${info.release_date || 'N/A'}</p>
            </div>
        `;

        // 尝试获取分析数据
        try {
            const analysisResponse = await fetch(`${API_BASE}/pdb/analyze/${pdbId}`);
            const analysis = await analysisResponse.json();

            if (!analysis.error) {
                html += `
                    <div class="detail-section">
                        <h3>🧪 物化性质</h3>
                        <p><strong>链数:</strong> ${analysis.num_chains || 'N/A'}</p>
                        <p><strong>残基数:</strong> ${analysis.num_residues || 'N/A'}</p>
                        <p><strong>原子数:</strong> ${analysis.num_atoms || 'N/A'}</p>
                    </div>
                `;

                if (analysis.secondary_structure) {
                    html += `
                        <div class="detail-section">
                            <h3>🔗 二级结构</h3>
                            <p><strong>螺旋:</strong> ${analysis.secondary_structure.helix || 'N/A'}</p>
                            <p><strong>β折叠:</strong> ${analysis.secondary_structure.beta_sheet || 'N/A'}</p>
                            <p><strong>线圈:</strong> ${analysis.secondary_structure.coil || 'N/A'}</p>
                        </div>
                    `;
                }
            }
        } catch (e) {
            console.log('分析数据获取失败:', e);
        }

        // 添加外部链接
        html += `
            <div class="detail-section">
                <h3>🔗 外部链接</h3>
                <p><a href="https://www.rcsb.org/structure/${pdbId}" target="_blank">RCSB PDB 官网 →</a></p>
                <p><a href="https://www.rcsb.org/3d-view/${pdbId}" target="_blank">RCSB 3D 查看器 →</a></p>
                <p><a href="https://molstar.org/viewer/?pdb-id=${pdbId}" target="_blank">Molstar 查看器 →</a></p>
            </div>
        `;

        content.innerHTML = html;

    } catch (error) {
        content.innerHTML = `<p>加载详情失败: ${error.message}</p>`;
    }
}

// 加载 3D 查看器
function load3DViewer(pdbId) {
    const container = document.getElementById('viewer3d');
    container.innerHTML = ''; // 清空

    // 创建新的查看器
    viewer = $3Dmol.createViewer(container, {
        backgroundColor: '#1a1a2e'
    });

    // 从 RCSB 加载 PDB 数据
    $3Dmol.download(`pdb:${pdbId}`, viewer, {}, function() {
        updateViewer();
        viewer.zoomTo();
        viewer.render();
    });
}

// 更新查看器样式
function updateViewer() {
    if (!viewer) return;

    const style = document.getElementById('styleSelect').value;
    const color = document.getElementById('colorSelect').value;

    viewer.setStyle({}, {}); // 清除所有样式

    let styleObj = {};
    let colorScheme = {};

    // 颜色方案
    switch (color) {
        case 'spectrum':
            colorScheme = { color: 'spectrum' };
            break;
        case 'chain':
            colorScheme = { colorscheme: 'chain' };
            break;
        case 'ss':
            colorScheme = { colorscheme: 'ssJmol' };
            break;
    }

    // 显示样式
    switch (style) {
        case 'cartoon':
            styleObj = { cartoon: colorScheme };
            break;
        case 'stick':
            styleObj = { stick: { ...colorScheme, radius: 0.2 } };
            break;
        case 'sphere':
            styleObj = { sphere: { ...colorScheme, radius: 0.5 } };
            break;
        case 'surface':
            styleObj = { cartoon: colorScheme };
            viewer.addSurface($3Dmol.VDW, { opacity: 0.7, color: 'white' });
            break;
    }

    viewer.setStyle({}, styleObj);
    viewer.render();
}

// 重置查看器
function resetViewer() {
    if (!viewer) return;
    viewer.zoomTo();
    viewer.render();
}

// 生成报告（基因搜索）
async function generateReport(geneName) {
    const content = document.getElementById('reportContent');
    content.textContent = '正在生成报告...';

    try {
        const response = await fetch(`${API_BASE}/report?gene_name=${encodeURIComponent(geneName)}`);
        const data = await response.json();

        if (data.error) {
            throw new Error(data.error);
        }

        currentReport = data.report;
        content.textContent = data.report;
    } catch (error) {
        content.textContent = `生成报告失败: ${error.message}`;
    }
}

// 生成报告（PDB 搜索）
async function generateReportForPdb(pdbId) {
    const content = document.getElementById('reportContent');
    content.textContent = '正在生成报告...';

    try {
        const response = await fetch(`${API_BASE}/report?pdb_ids=${pdbId}`);
        const data = await response.json();

        if (data.error) {
            throw new Error(data.error);
        }

        currentReport = data.report;
        content.textContent = data.report;
    } catch (error) {
        content.textContent = `生成报告失败: ${error.message}`;
    }
}

// 下载报告
function downloadReport() {
    if (!currentReport) {
        showError('没有可下载的报告');
        return;
    }

    const blob = new Blob([currentReport], { type: 'text/markdown' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `pdb_analysis_report_${new Date().toISOString().slice(0, 10)}.md`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

// UI 辅助函数
function showLoading(show) {
    document.getElementById('loading').classList.toggle('hidden', !show);
}

function showResults() {
    document.getElementById('results').classList.remove('hidden');
}

function hideResults() {
    document.getElementById('results').classList.add('hidden');
}

function showError(message) {
    const errorMsg = document.getElementById('errorMsg');
    document.getElementById('errorText').textContent = message;
    errorMsg.classList.remove('hidden');

    // 5秒后自动隐藏
    setTimeout(hideError, 5000);
}

function hideError() {
    document.getElementById('errorMsg').classList.add('hidden');
}

