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

        // 加载高级分析
        loadAdvancedAnalysis(pdbId);

        // 加载序列分析
        loadSequenceAnalysis(pdbId);

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
    content.innerHTML = '正在生成报告...';

    try {
        const response = await fetch(`${API_BASE}/report?gene_name=${encodeURIComponent(geneName)}`);
        const data = await response.json();

        if (data.error) {
            throw new Error(data.error);
        }

        currentReport = data.report;
        // 使用 marked 解析 Markdown
        content.innerHTML = marked.parse(data.report);
    } catch (error) {
        content.innerHTML = `生成报告失败: ${error.message}`;
    }
}

// 生成报告（PDB 搜索）
async function generateReportForPdb(pdbId) {
    const content = document.getElementById('reportContent');
    content.innerHTML = '正在生成报告...';

    try {
        const response = await fetch(`${API_BASE}/report?pdb_ids=${pdbId}`);
        const data = await response.json();

        if (data.error) {
            throw new Error(data.error);
        }

        currentReport = data.report;
        // 使用 marked 解析 Markdown
        content.innerHTML = marked.parse(data.report);
    } catch (error) {
        content.innerHTML = `生成报告失败: ${error.message}`;
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

// ==================== 高级结构分析 ====================
async function loadAdvancedAnalysis(pdbId) {
    const container = document.getElementById('advancedAnalysis');
    const content = document.getElementById('advancedContent');

    container.classList.remove('hidden');
    content.innerHTML = '<p>正在加载高级分析数据...</p>';

    try {
        const response = await fetch(`${API_BASE}/pdb/analyze-advanced/${pdbId}`);
        const data = await response.json();

        if (data.error) {
            content.innerHTML = `<p>加载失败: ${data.error}</p>`;
            return;
        }

        let html = '';

        // 二硫键
        html += `
            <div class="analysis-section">
                <h3>🔗 二硫键</h3>
                <p><strong>数量:</strong> ${data.disulfide_bonds?.count || 0}</p>
        `;
        if (data.disulfide_bonds?.bonds?.length > 0) {
            html += '<ul>';
            data.disulfide_bonds.bonds.forEach(bond => {
                html += `<li>${bond.cys1} ↔ ${bond.cys2} (${bond.distance}Å)</li>`;
            });
            html += '</ul>';
        }
        html += '</div>';

        // 盐桥
        html += `
            <div class="analysis-section">
                <h3>⚡ 盐桥</h3>
                <p><strong>数量:</strong> ${data.salt_bridges?.count || 0}</p>
        `;
        if (data.salt_bridges?.bridges?.length > 0) {
            html += '<ul class="salt-bridge-list">';
            data.salt_bridges.bridges.slice(0, 10).forEach(bridge => {
                html += `<li>${bridge.positive} ↔ ${bridge.negative} (${bridge.distance}Å)</li>`;
            });
            if (data.salt_bridges.bridges.length > 10) {
                html += `<li>...及其他 ${data.salt_bridges.bridges.length - 10} 个</li>`;
            }
            html += '</ul>';
        }
        html += '</div>';

        // 氢键
        html += `
            <div class="analysis-section">
                <h3>💧 氢键</h3>
                <p><strong>主链氢键数:</strong> ${data.hydrogen_bonds?.backbone_hbonds || 'N/A'}</p>
            </div>
        `;

        // SASA
        if (data.sasa_per_chain && !data.sasa_per_chain.error) {
            html += `
                <div class="analysis-section">
                    <h3>🌊 溶剂可及表面积 (SASA)</h3>
                    <table class="data-table">
                        <tr><th>链</th><th>SASA (Å²)</th></tr>
            `;
            for (const [chain, sasa] of Object.entries(data.sasa_per_chain)) {
                html += `<tr><td>链 ${chain}</td><td>${sasa}</td></tr>`;
            }
            html += '</table></div>';
        }

        // 疏水/亲水比例
        if (data.hydrophobicity_per_chain) {
            html += `
                <div class="analysis-section">
                    <h3>💦 疏水/亲水残基比例</h3>
                    <table class="data-table">
                        <tr><th>链</th><th>疏水残基</th><th>亲水残基</th><th>比例</th></tr>
            `;
            for (const [chain, info] of Object.entries(data.hydrophobicity_per_chain)) {
                html += `
                    <tr>
                        <td>链 ${chain}</td>
                        <td>${info.hydrophobic_count} (${info.hydrophobic_ratio}%)</td>
                        <td>${info.hydrophilic_count} (${info.hydrophilic_ratio}%)</td>
                        <td>
                            <div class="ratio-bar">
                                <div class="hydrophobic" style="width: ${info.hydrophobic_ratio}%"></div>
                                <div class="hydrophilic" style="width: ${info.hydrophilic_ratio}%"></div>
                            </div>
                        </td>
                    </tr>
                `;
            }
            html += '</table></div>';
        }

        content.innerHTML = html;

    } catch (error) {
        content.innerHTML = `<p>加载高级分析失败: ${error.message}</p>`;
    }
}

// ==================== 突变影响分析 ====================
async function analyzeMutation() {
    const mutationInput = document.getElementById('mutationInput').value.trim();
    const resultContainer = document.getElementById('mutationResult');

    if (!currentPdbId) {
        showError('请先选择一个 PDB 结构');
        return;
    }

    if (!mutationInput) {
        showError('请输入突变信息，格式: A:K33E');
        return;
    }

    resultContainer.classList.remove('hidden');
    resultContainer.innerHTML = '<p>正在分析突变影响...</p>';

    try {
        const response = await fetch(
            `${API_BASE}/pdb/mutation?pdb_id=${currentPdbId}&mutation=${encodeURIComponent(mutationInput)}`
        );
        const data = await response.json();

        if (data.error) {
            resultContainer.innerHTML = `<p class="error">${data.error}</p>`;
            return;
        }

        const impact = data.impact_assessment;
        const impactClass = impact.level === '高' ? 'high' : (impact.level === '中' ? 'medium' : 'low');

        let html = `
            <div class="mutation-summary">
                <h3>突变: ${data.mutation}</h3>
                <div class="impact-badge ${impactClass}">影响程度: ${impact.level}</div>
                <p>${impact.description}</p>
            </div>
            
            <div class="mutation-details">
                <div class="aa-comparison">
                    <div class="aa-box wt">
                        <h4>野生型 (${data.wild_type.aa})</h4>
                        <p><strong>名称:</strong> ${data.wild_type.name}</p>
                        <p><strong>电荷:</strong> ${data.wild_type.charge}</p>
                        <p><strong>体积:</strong> ${data.wild_type.volume}Å³</p>
                        <p><strong>疏水性:</strong> ${data.wild_type.hydrophobic ? '是' : '否'}</p>
                    </div>
                    <div class="aa-arrow">→</div>
                    <div class="aa-box mut">
                        <h4>突变型 (${data.mutant.aa})</h4>
                        <p><strong>名称:</strong> ${data.mutant.name}</p>
                        <p><strong>电荷:</strong> ${data.mutant.charge}</p>
                        <p><strong>体积:</strong> ${data.mutant.volume}Å³</p>
                        <p><strong>疏水性:</strong> ${data.mutant.hydrophobic ? '是' : '否'}</p>
                    </div>
                </div>
                
                <div class="changes-summary">
                    <h4>变化摘要</h4>
                    <ul>
                        <li>电荷变化: ${data.changes.charge_change > 0 ? '+' : ''}${data.changes.charge_change}</li>
                        <li>体积变化: ${data.changes.volume_change > 0 ? '+' : ''}${data.changes.volume_change}Å³</li>
                        <li>疏水性变化: ${data.changes.hydrophobicity_change ? '是' : '否'}</li>
                        <li>极性变化: ${data.changes.polarity_change ? '是' : '否'}</li>
                    </ul>
                </div>
        `;

        if (impact.reasons && impact.reasons.length > 0) {
            html += `
                <div class="impact-reasons">
                    <h4>影响原因</h4>
                    <ul>
                        ${impact.reasons.map(r => `<li>${r}</li>`).join('')}
                    </ul>
                </div>
            `;
        }

        if (data.structural_context) {
            html += `
                <div class="structural-context">
                    <h4>结构上下文</h4>
                    <p><strong>结构中该位置残基:</strong> ${data.structural_context.found_residue || 'N/A'}</p>
                    ${data.structural_context.secondary_structure ? `<p><strong>二级结构:</strong> ${data.structural_context.secondary_structure}</p>` : ''}
                    ${data.structural_context.warning ? `<p class="warning">⚠️ ${data.structural_context.warning}</p>` : ''}
                </div>
            `;
        }

        html += '</div>';
        resultContainer.innerHTML = html;

    } catch (error) {
        resultContainer.innerHTML = `<p class="error">分析失败: ${error.message}</p>`;
    }
}

// ==================== 序列组成分析 ====================
let sequenceCharts = {};

async function loadSequenceAnalysis(pdbId) {
    const container = document.getElementById('sequenceAnalysis');
    const content = document.getElementById('sequenceContent');

    container.classList.remove('hidden');
    content.innerHTML = '<p>正在加载序列分析数据...</p>';

    try {
        const response = await fetch(`${API_BASE}/pdb/sequence-composition/${pdbId}`);
        const data = await response.json();

        if (data.error) {
            content.innerHTML = `<p>加载失败: ${data.error}</p>`;
            return;
        }

        let html = '';

        // 为每条链创建图表
        for (const [chainId, chainData] of Object.entries(data.chains || {})) {
            html += `
                <div class="chain-analysis">
                    <h3>链 ${chainId} (${chainData.length} 个残基)</h3>
                    <div class="category-stats">
                        <span class="stat-item positive">正电荷: ${chainData.category_statistics.charged_positive_pct}%</span>
                        <span class="stat-item negative">负电荷: ${chainData.category_statistics.charged_negative_pct}%</span>
                        <span class="stat-item hydrophobic">疏水: ${chainData.category_statistics.hydrophobic_pct}%</span>
                        <span class="stat-item polar">极性: ${chainData.category_statistics.polar_uncharged_pct}%</span>
                        <span class="stat-item aromatic">芳香: ${chainData.category_statistics.aromatic_pct}%</span>
                    </div>
                    <canvas id="chart-${pdbId}-${chainId}" height="200"></canvas>
                </div>
            `;
        }

        content.innerHTML = html;

        // 绘制图表
        for (const [chainId, chainData] of Object.entries(data.chains || {})) {
            createAminoAcidChart(`chart-${pdbId}-${chainId}`, chainData);
        }

    } catch (error) {
        content.innerHTML = `<p>加载序列分析失败: ${error.message}</p>`;
    }
}

function createAminoAcidChart(canvasId, chainData) {
    const ctx = document.getElementById(canvasId);
    if (!ctx) return;

    // 销毁之前的图表
    if (sequenceCharts[canvasId]) {
        sequenceCharts[canvasId].destroy();
    }

    const aminoAcids = Object.keys(chainData.amino_acid_percentages);
    const percentages = Object.values(chainData.amino_acid_percentages);

    // 根据氨基酸属性设置颜色
    const colors = aminoAcids.map(aa => {
        if (['K', 'R', 'H'].includes(aa)) return '#3498db'; // 正电荷 - 蓝色
        if (['D', 'E'].includes(aa)) return '#e74c3c'; // 负电荷 - 红色
        if (['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'].includes(aa)) return '#f39c12'; // 疏水 - 橙色
        if (['S', 'T', 'N', 'Q', 'Y', 'C'].includes(aa)) return '#2ecc71'; // 极性 - 绿色
        return '#9b59b6'; // 其他 - 紫色
    });

    sequenceCharts[canvasId] = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: aminoAcids,
            datasets: [{
                label: '氨基酸占比 (%)',
                data: percentages,
                backgroundColor: colors,
                borderColor: colors.map(c => c),
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            plugins: {
                legend: {
                    display: false
                },
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            return `${context.parsed.y.toFixed(2)}%`;
                        }
                    }
                }
            },
            scales: {
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: '百分比 (%)'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: '氨基酸'
                    }
                }
            }
        }
    });
}

// ==================== UniProt 序列比对 ====================
async function alignWithUniprot() {
    const uniprotId = document.getElementById('uniprotIdInput').value.trim();
    const container = document.getElementById('alignmentAnalysis');
    const result = document.getElementById('alignmentResult');

    if (!currentPdbId) {
        showError('请先选择一个 PDB 结构');
        return;
    }

    container.classList.remove('hidden');
    result.innerHTML = '<p>正在进行序列比对...</p>';

    try {
        let url = `${API_BASE}/pdb/align-uniprot/${currentPdbId}`;
        if (uniprotId) {
            url += `?uniprot_id=${encodeURIComponent(uniprotId)}`;
        }

        const response = await fetch(url);
        const data = await response.json();

        if (data.error) {
            result.innerHTML = `<p class="error">${data.error}</p>`;
            return;
        }

        let html = `
            <div class="alignment-summary">
                <h3>比对结果</h3>
                <p><strong>PDB ID:</strong> ${data.pdb_id}</p>
                <p><strong>UniProt ID:</strong> ${data.uniprot_id}</p>
                <p><strong>UniProt 序列长度:</strong> ${data.uniprot_length} 残基</p>
            </div>
            
            <div class="chain-alignments">
        `;

        for (const [chainId, alignment] of Object.entries(data.chain_alignments || {})) {
            html += `
                <div class="chain-alignment">
                    <h4>链 ${chainId}</h4>
                    <div class="alignment-stats">
                        <div class="stat">
                            <span class="label">PDB 长度</span>
                            <span class="value">${alignment.pdb_length}</span>
                        </div>
                        <div class="stat">
                            <span class="label">序列一致性</span>
                            <span class="value identity">${alignment.identity_percent}%</span>
                        </div>
                        <div class="stat">
                            <span class="label">覆盖率</span>
                            <span class="value">${alignment.coverage_percent}%</span>
                        </div>
                    </div>
            `;

            if (alignment.missing_regions && alignment.missing_regions.length > 0) {
                html += `
                    <div class="missing-regions">
                        <h5>缺失区段</h5>
                        <ul>
                            ${alignment.missing_regions.map(r => 
                                `<li>位置 ${r.start}-${r.end} (${r.length} 残基)</li>`
                            ).join('')}
                        </ul>
                    </div>
                `;
            } else {
                html += '<p class="no-missing">✅ 无缺失区段</p>';
            }

            html += '</div>';
        }

        html += '</div>';
        result.innerHTML = html;

    } catch (error) {
        result.innerHTML = `<p class="error">比对失败: ${error.message}</p>`;
    }
}
