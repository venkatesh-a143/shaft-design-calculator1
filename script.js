// Toggle functions for UI
function toggleProblemType() {
    const problemType = document.querySelector('input[name="problemType"]:checked').value;
    document.getElementById('udlSection').style.display = problemType === 'udl' ? 'block' : 'none';
}

function toggleStressMethod() {
    const method = document.querySelector('input[name="stressMethod"]:checked').value;
    document.getElementById('directStressInputs').style.display = method === 'direct' ? 'block' : 'none';
    document.getElementById('materialInputs').style.display = method === 'material' ? 'block' : 'none';
}

function toggleCmCt() {
    const useCmCt = document.getElementById('useCmCt').checked;
    document.getElementById('cmCtInputs').style.display = useCmCt ? 'block' : 'none';
}

function toggleShaftType() {
    const shaftType = document.querySelector('input[name="shaftType"]:checked').value;
    document.getElementById('hollowInputs').style.display = shaftType === 'hollow' ? 'block' : 'none';
}

// Evaluate fraction input
function evalFraction(str) {
    if (str.includes('/')) {
        const parts = str.split('/');
        return parseFloat(parts[0]) / parseFloat(parts[1]);
    }
    return parseFloat(str);
}

// Main calculation function
function calculateShaft() {
    // Get inputs
    const problemType = document.querySelector('input[name="problemType"]:checked').value;
    const P = parseFloat(document.getElementById('power').value);
    const N = parseFloat(document.getElementById('speed').value);
    
    // Calculate torque (N-mm)
    const T = (9.55e6 * P) / N;
    
    // Get stress parameters
    const stressMethod = document.querySelector('input[name="stressMethod"]:checked').value;
    let sigmaMax, tauMax;
    
    if (stressMethod === 'direct') {
        tauMax = parseFloat(document.getElementById('tauAllow').value);
        sigmaMax = 2 * tauMax;
    } else {
        const sigmaUt = parseFloat(document.getElementById('sigmaUt').value);
        const sigmaYt = parseFloat(document.getElementById('sigmaYt').value);
        const nPrime = parseFloat(document.getElementById('nPrime').value);
        
        const sigmaU = sigmaUt / nPrime;
        const sigmaY = sigmaYt / nPrime;
        
        const sigmaMax1 = 0.36 * sigmaU;
        const sigmaMax2 = 0.6 * sigmaY;
        sigmaMax = Math.min(sigmaMax1, sigmaMax2);
        
        const tauMax1 = 0.18 * sigmaU;
        const tauMax2 = 0.3 * sigmaY;
        tauMax = Math.min(tauMax1, tauMax2);
    }
    
    // Get Cm and Ct
    const useCmCt = document.getElementById('useCmCt').checked;
    const Cm = useCmCt ? parseFloat(document.getElementById('cm').value) : 1.5;
    const Ct = useCmCt ? parseFloat(document.getElementById('ct').value) : 1.0;
    
    // Calculate shaft configuration
    let L, bearing1, bearing2, maxBM;
    let configData = {};
    
    if (problemType === 'udl') {
        L = parseFloat(document.getElementById('totalLength').value);
        const lUdl = parseFloat(document.getElementById('udlLength').value);
        const w = parseFloat(document.getElementById('udlIntensity').value);
        
        bearing1 = (L - lUdl) / 2;
        bearing2 = bearing1 + lUdl;
        
        const WTotal = w * lUdl;
        const RA = WTotal / 2;
        const RB = WTotal / 2;
        
        // Maximum bending moment for UDL (at center)
        maxBM = (w * lUdl * lUdl) / 8;
        
        configData = {
            type: 'udl',
            L: L,
            bearing1: bearing1,
            bearing2: bearing2,
            w: w,
            lUdl: lUdl,
            WTotal: WTotal,
            RA: RA,
            RB: RB
        };
    } else {
        // Point load
        L = 600;
        bearing1 = 100;
        bearing2 = 500;
        const gearPos = 300;
        
        const Dg = 200;
        const Ft = (2 * T) / Dg;
        const Fr = Ft * Math.tan(20 * Math.PI / 180);
        
        const RBh = (Ft * (gearPos - bearing1)) / (bearing2 - bearing1);
        const RAh = Ft - RBh;
        
        const RBv = (Fr * (gearPos - bearing1)) / (bearing2 - bearing1);
        const RAv = Fr - RBv;
        
        const BMh = RAh * (gearPos - bearing1);
        const BMv = RAv * (gearPos - bearing1);
        
        maxBM = Math.sqrt(BMh * BMh + BMv * BMv);
        
        configData = {
            type: 'point',
            L: L,
            bearing1: bearing1,
            bearing2: bearing2,
            gearPos: gearPos,
            Ft: Ft,
            Fr: Fr,
            RAh: RAh,
            RBh: RBh,
            RAv: RAv,
            RBv: RBv
        };
    }
    
    // Calculate shaft diameter
    const shaftType = document.querySelector('input[name="shaftType"]:checked').value;
    let diameter, outerDia, innerDia;
    
    const term1 = Cm * maxBM;
    const term2 = Math.sqrt(Math.pow(Cm * maxBM, 2) + Math.pow(Ct * T, 2));
    
    if (shaftType === 'solid') {
        // Method 1: Normal stress
        const dNormal = Math.pow((16 / (Math.PI * sigmaMax)) * (term1 + term2), 1/3);
        
        // Method 2: Shear stress
        const dShear = Math.pow((16 / (Math.PI * tauMax)) * term2, 1/3);
        
        diameter = Math.max(dNormal, dShear);
        
        // Standard size
        const standardSizes = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100];
        diameter = standardSizes.find(s => s >= diameter) || Math.ceil(diameter / 10) * 10;
        
        outerDia = diameter;
        innerDia = 0;
    } else {
        // Hollow shaft
        const kStr = document.getElementById('diameterRatio').value;
        const k = evalFraction(kStr);
        
        // Method 1: Normal stress
        const doNormal = Math.pow((16 / (Math.PI * sigmaMax)) * (term1 + term2) * (1 / (1 - Math.pow(k, 4))), 1/3);
        
        // Method 2: Shear stress
        const doShear = Math.pow((16 / (Math.PI * tauMax)) * term2 * (1 / (1 - Math.pow(k, 4))), 1/3);
        
        outerDia = Math.max(doNormal, doShear);
        
        // Standard size
        const standardSizes = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100];
        outerDia = standardSizes.find(s => s >= outerDia) || Math.ceil(outerDia / 10) * 10;
        innerDia = k * outerDia;
        
        diameter = outerDia;
    }
    
    // Display results
    displayResults(T, maxBM, diameter, outerDia, innerDia, shaftType, Cm, Ct, sigmaMax, tauMax);
    
    // Generate diagrams
    generateDiagrams(configData, T, maxBM);
}

function displayResults(T, maxBM, diameter, outerDia, innerDia, shaftType, Cm, Ct, sigmaMax, tauMax) {
    document.getElementById('results').style.display = 'block';
    document.getElementById('diagrams').style.display = 'block';
    
    document.getElementById('torqueResult').textContent = `${T.toFixed(2)} N-mm (${(T/1000).toFixed(2)} N-m)`;
    document.getElementById('momentResult').textContent = `${maxBM.toFixed(2)} N-mm (${(maxBM/1000).toFixed(2)} N-m)`;
    
    if (shaftType === 'solid') {
        document.getElementById('diameterResult').textContent = `d = ${diameter.toFixed(2)} mm`;
    } else {
        document.getElementById('diameterResult').textContent = `do = ${outerDia.toFixed(2)} mm, di = ${innerDia.toFixed(2)} mm`;
    }
    
    // Detailed results
    const detailedHTML = `
        <h3>Calculation Details</h3>
        <p><strong>Combined Shock Factors:</strong> Cm = ${Cm.toFixed(2)}, Ct = ${Ct.toFixed(2)}</p>
        <p><strong>Maximum Normal Stress (σ_max):</strong> ${sigmaMax.toFixed(2)} MPa</p>
        <p><strong>Maximum Shear Stress (τ_max):</strong> ${tauMax.toFixed(2)} MPa</p>
        <p><strong>Transmitted Torque (T):</strong> ${T.toFixed(2)} N-mm</p>
        <p><strong>Maximum Bending Moment (M):</strong> ${maxBM.toFixed(2)} N-mm</p>
        ${shaftType === 'solid' ? 
            `<p><strong>✓ Solid Shaft Diameter:</strong> ${diameter.toFixed(2)} mm</p>` :
            `<p><strong>✓ Hollow Shaft:</strong> Outer Diameter = ${outerDia.toFixed(2)} mm, Inner Diameter = ${innerDia.toFixed(2)} mm</p>`
        }
    `;
    
    document.getElementById('detailedResults').innerHTML = detailedHTML;
    
    // Scroll to results
    document.getElementById('results').scrollIntoView({ behavior: 'smooth' });
}

function generateDiagrams(configData, T, maxBM) {
    const n = 500;
    const x = Array.from({length: n}, (_, i) => (i / (n - 1)) * configData.L);
    
    // Generate SF and BM data for HORIZONTAL plane
    const SF_h = [];
    const BM_h = [];
    
    // Generate SF and BM data for VERTICAL plane
    const SF_v = [];
    const BM_v = [];
    
    if (configData.type === 'udl') {
        // UDL - Same for both horizontal and vertical planes
        for (let i = 0; i < n; i++) {
            if (x[i] < configData.bearing1 || x[i] > configData.bearing2) {
                SF_h.push(0);
                BM_h.push(0);
                SF_v.push(0);
                BM_v.push(0);
            } else {
                const dist = x[i] - configData.bearing1;
                const sf = configData.RA - configData.w * dist;
                const bm = configData.RA * dist - configData.w * dist * dist / 2;
                
                SF_h.push(sf);
                BM_h.push(bm);
                SF_v.push(sf);
                BM_v.push(bm);
            }
        }
    } else {
        // Point load - separate horizontal and vertical
        const gearPos = configData.gearPos;
        
        for (let i = 0; i < n; i++) {
            // Horizontal plane
            if (x[i] < configData.bearing1) {
                SF_h.push(0);
                BM_h.push(0);
            } else if (x[i] >= configData.bearing1 && x[i] < gearPos) {
                SF_h.push(configData.RAh);
                BM_h.push(configData.RAh * (x[i] - configData.bearing1));
            } else if (x[i] >= gearPos && x[i] <= configData.bearing2) {
                SF_h.push(configData.RAh - configData.Ft);
                BM_h.push(configData.RAh * (x[i] - configData.bearing1) - configData.Ft * (x[i] - gearPos));
            } else {
                SF_h.push(0);
                BM_h.push(0);
            }
            
            // Vertical plane
            if (x[i] < configData.bearing1) {
                SF_v.push(0);
                BM_v.push(0);
            } else if (x[i] >= configData.bearing1 && x[i] < gearPos) {
                SF_v.push(configData.RAv);
                BM_v.push(configData.RAv * (x[i] - configData.bearing1));
            } else if (x[i] >= gearPos && x[i] <= configData.bearing2) {
                SF_v.push(configData.RAv - configData.Fr);
                BM_v.push(configData.RAv * (x[i] - configData.bearing1) - configData.Fr * (x[i] - gearPos));
            } else {
                SF_v.push(0);
                BM_v.push(0);
            }
        }
    }
    
    // Calculate resultant bending moment
    const BM_resultant = [];
    for (let i = 0; i < n; i++) {
        BM_resultant.push(Math.sqrt(BM_h[i] * BM_h[i] + BM_v[i] * BM_v[i]));
    }
    
    // Generate torque data
    const torqueData = configData.type === 'udl' ? 
        Array(n).fill(T) : 
        x.map(xi => xi >= configData.gearPos ? T : 0);
    
    // Clear previous diagrams
    document.getElementById('diagrams').innerHTML = `
        <h2>📈 Horizontal Loading Diagram (HLD)</h2>
        <div id="sfdDiagram_h"></div>
        <div id="bmdDiagram_h"></div>
        
        <h2>📈 Vertical Loading Diagram (VLD)</h2>
        <div id="sfdDiagram_v"></div>
        <div id="bmdDiagram_v"></div>
        
        <h2>📈 Resultant Bending Moment</h2>
        <div id="bmdDiagram_resultant"></div>
        
        <h2>📈 Torque Diagram</h2>
        <div id="torqueDiagram"></div>
    `;
    
    // Plot HLD - SFD
    Plotly.newPlot('sfdDiagram_h', [{
        x: x,
        y: SF_h,
        type: 'scatter',
        fill: 'tozeroy',
        name: 'Shear Force (Horizontal)',
        line: {color: 'rgb(55, 128, 191)', width: 3}
    }], {
        title: 'Shear Force Diagram - Horizontal Plane (SFD)',
        xaxis: {title: 'Distance (mm)'},
        yaxis: {title: 'Shear Force (N)'},
        showlegend: false
    });
    
    // Plot HLD - BMD
    const maxBM_h = Math.max(...BM_h.map(Math.abs));
    const maxBM_h_idx = BM_h.findIndex(bm => Math.abs(bm) === maxBM_h);
    
    Plotly.newPlot('bmdDiagram_h', [{
        x: x,
        y: BM_h,
        type: 'scatter',
        fill: 'tozeroy',
        name: 'Bending Moment (Horizontal)',
        line: {color: 'rgb(219, 64, 82)', width: 3}
    }], {
        title: 'Bending Moment Diagram - Horizontal Plane (BMD)',
        xaxis: {title: 'Distance (mm)'},
        yaxis: {title: 'Bending Moment (N-mm)'},
        showlegend: false,
        annotations: [{
            x: x[maxBM_h_idx],
            y: BM_h[maxBM_h_idx],
            text: `M_max = ${BM_h[maxBM_h_idx].toFixed(0)} N-mm`,
            showarrow: true,
            arrowhead: 2,
            ax: 0,
            ay: -40
        }]
    });
    
    // Plot VLD - SFD
    Plotly.newPlot('sfdDiagram_v', [{
        x: x,
        y: SF_v,
        type: 'scatter',
        fill: 'tozeroy',
        name: 'Shear Force (Vertical)',
        line: {color: 'rgb(55, 128, 191)', width: 3}
    }], {
        title: 'Shear Force Diagram - Vertical Plane (SFD)',
        xaxis: {title: 'Distance (mm)'},
        yaxis: {title: 'Shear Force (N)'},
        showlegend: false
    });
    
    // Plot VLD - BMD
    const maxBM_v = Math.max(...BM_v.map(Math.abs));
    const maxBM_v_idx = BM_v.findIndex(bm => Math.abs(bm) === maxBM_v);
    
    Plotly.newPlot('bmdDiagram_v', [{
        x: x,
        y: BM_v,
        type: 'scatter',
        fill: 'tozeroy',
        name: 'Bending Moment (Vertical)',
        line: {color: 'rgb(219, 64, 82)', width: 3}
    }], {
        title: 'Bending Moment Diagram - Vertical Plane (BMD)',
        xaxis: {title: 'Distance (mm)'},
        yaxis: {title: 'Bending Moment (N-mm)'},
        showlegend: false,
        annotations: [{
            x: x[maxBM_v_idx],
            y: BM_v[maxBM_v_idx],
            text: `M_max = ${BM_v[maxBM_v_idx].toFixed(0)} N-mm`,
            showarrow: true,
            arrowhead: 2,
            ax: 0,
            ay: -40
        }]
    });
    
    // Plot Resultant BMD
    const maxBM_res = Math.max(...BM_resultant);
    const maxBM_res_idx = BM_resultant.indexOf(maxBM_res);
    
    Plotly.newPlot('bmdDiagram_resultant', [{
        x: x,
        y: BM_resultant,
        type: 'scatter',
        fill: 'tozeroy',
        name: 'Resultant Bending Moment',
        line: {color: 'rgb(148, 0, 211)', width: 3}
    }], {
        title: 'Resultant Bending Moment Diagram (√(M_h² + M_v²))',
        xaxis: {title: 'Distance (mm)'},
        yaxis: {title: 'Bending Moment (N-mm)'},
        showlegend: false,
        annotations: [{
            x: x[maxBM_res_idx],
            y: BM_resultant[maxBM_res_idx],
            text: `M_max = ${maxBM_res.toFixed(0)} N-mm`,
            showarrow: true,
            arrowhead: 2,
            ax: 0,
            ay: -40
        }]
    });
    
    // Plot Torque Diagram
    Plotly.newPlot('torqueDiagram', [{
        x: x,
        y: torqueData,
        type: 'scatter',
        fill: 'tozeroy',
        name: 'Torque',
        line: {color: 'rgb(50, 171, 96)', width: 3}
    }], {
        title: 'Torque Diagram',
        xaxis: {title: 'Distance (mm)'},
        yaxis: {title: 'Torque (N-mm)'},
        showlegend: false
    });
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', function() {
    console.log('Shaft Design Calculator loaded');
});
