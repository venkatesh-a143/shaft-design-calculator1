// Toggle functions for UI
function toggleProblemType() {
    const val = document.querySelector('input[name="problemType"]:checked').value;
    document.getElementById('udlSection').style.display = val === 'udl' ? 'block' : 'none';
    document.getElementById('twoGearSection').style.display = val === 'twogear' ? 'block' : 'none';
    document.getElementById('pulleyGearSection').style.display = val === 'pulleygear' ? 'block' : 'none';
}

function toggleStressMethod() {
    const method = document.querySelector('input[name="stressMethod"]:checked').>Allowable Shear Stress: τ = ${tauMax.toFixed(2)} MPa</p>`;
        details += `<p>Allowable Normal Stress: σ = ${sigmaMax.toFixed(2)} MPa</p>`;
    } else {
        const sigmaUt = parseFloat(document.getElementById('sigmaUt').value);
        const sigmaYt = parseFloat(document.getElementById('sigmaYt').value);
        const nPrime = parseFloat(document.getElementById('nPrime').value);

        const sigmaU = sigmaUt / nPrime;
        const sigmaY = sigmaYt / nPrime;

        details += `<p>σ_u = σ_ut/n' = ${sigmaUt}/${nPrime} = ${sigmaU.toFixed(2)} MPa</p>`;
        details += `<p>σ_y = σ_yt/n' = ${sigmaYt}/${nPrime} = ${sigmaY.toFixed(2)} MPa</p>`;

        const sm1 = 0.36 * sigmaU;
        const sm2 = 0.6 * sigmaY;
        sigmaMax = Math.min(sm1, sm2);

        details += `<p>σ_max = 0.36×σ_u = ${sm1.toFixed(2)} MPa</p>`;
        details += `<p>σ_max = 0.6×σ_y = ${sm2.toFixed(2)} MPa</p>`;
        details += `<p><strong>σ_max = ${sigmaMax.toFixed(2)} MPa (least value)</strong></p>`;

        const tm1 = 0.18value;
    document.getElementById('directStressInputs').style.display = method === 'direct' ? 'block' : 'none';
    document.getElementById('materialInputs').style.display = method === 'material' ? 'block' : 'none';
    document.getElementById('workingStressInputs').style.display = method === 'working' ? 'block' : 'none';
}

function toggleCmCt() {
    document.getElementById('cmCtInputs').style.display = document.getElementById('useCmCt').checked ? 'block' : 'none';
}

function toggleShaftType() {
    document.getElementById('hollowInputs').style.display = document.querySelector('input[name="shaftType"]:checked').value === 'hollow' ? 'block' : 'none';
}

function evalFraction(str) {
    if (str.includes('/')) {
        const parts = str.split('/');
        return parseFloat(parts[0]) / parseFloat(parts[1]);
    }
    return parseFloat(str);
}

// Main calculation function
function calculateShaft() {
    const problemType = document.querySelector('input[name="problemType"]:checked').value;
    const P = parseFloat(document.getElementById('power').value);
    const N = parseFloat(document.getElementById('speed').value);
    const T = (9.55e6 * P) / N;

    // Get stress parameters
    const stressMethod = document.querySelector('input[name="stressMethod"]:checked').value;
    let sigmaMax, tauMax;

    if (stressMethod === 'direct') {
        tauMax = parse * sigmaU;
        const tm2 = 0.3 * sigmaY;
        tauMax = Math.min(tm1, tm2);

        details += `<p>τ_max = 0.18×σ_u = ${tm1.toFixed(2)} MPa</p>`;
        details += `<p>τ_max = 0.3×σ_y = ${tm2.toFixed(2)} MPa</p>`;
        details += `<p><strong>τ_max = ${tauMax.toFixed(2)} MPa (least value)</strong></p>`;
    }

    if (hasKeyway) {
        const smOld = sigmaMax;
        const tmOld = tauMax;
        sigmaMax = 0.75 * sigmaMax;
        tauMax = 0.75 * tauMax;
        details += `<p><strong>With keyway (25% reduction):</strong></p>`;
        details += `<p>σ_max = 0.75 × ${smOld.toFixed(2)} = ${sigmaMax.toFixed(2)} MPa</p>`;
        details += `<p>τ_max = 0.75 × ${tmOld.toFixed(2)} = ${tauMax.toFixed(2)} MPa</p>`;
    }

    return { sigmaMax, tauMax, details };
}

// Main calculation
function calculateShaft() {
    const problemType = document.querySelector('input[name="problemType"]:checked').value;
    const P = parseFloat(document.getElementById('power').value);
    const N = parseFloat(document.getElementById('speed').value);
    const T = (9.55e6 * P) / N;

    const useCmCt = document.getElementByIdFloat(document.getElementById('tauAllow').value);
        sigmaMax = 2 * tauMax;
    } else if (stressMethod === 'material') {
        const sigmaUt = parseFloat(document.getElementById('sigmaUt').value);
        const sigmaYt = parseFloat(document.getElementById('sigmaYt').value);
        const nPrime = parseFloat(document.getElementById('nPrime').value);
        const sigmaU = sigmaUt / nPrime;
        const sigmaY = sigmaYt / nPrime;
        sigmaMax = Math.min(0.36 * sigmaU, 0.6 * sigmaY);
        tauMax = Math.min(0.18 * sigmaU, 0.3 * sigmaY);
    } else {
        sigmaMax = parseFloat(document.getElementById('sigmaWorking').value);
        tauMax = parseFloat(document.getElementById('tauWorking').value);
    }

    // Keyway reduction
    if (document.getElementById('hasKeyway').checked) {
        sigmaMax = 0.75 * sigmaMax;
        tauMax = 0.75 * tauMax;
    }

    // Get Cm and Ct
    const useCmCt = document.getElementById('useCmCt').checked;
    const Cm = useCmCt ? parseFloat(document.getElementById('cm').value) : 1.5;
    const Ct = useCmCt ? parseFloat(document.getElementById('ct').value) : 1.0;

    let result;

    if (problemType === 'point') {
        result = calcPointLoad(T);
    } else if (problemType === 'udl') {('useCmCt').checked;
    const Cm = useCmCt ? parseFloat(document.getElementById('cm').value) : 1.5;
    const Ct = useCmCt ? parseFloat(document.getElementById('ct').value) : 1.0;

    const stress = getStressValues();
    let configData;

    switch (problemType) {
        case 'point':
            configData = calcPointLoad(T);
            break;
        case 'udl':
            configData = calcUDL(T);
            break;
        case 'twogear':
            configData = calcTwoGear(T);
            break;
        case 'pulleygear':
            configData = calcPulleyGear(T);
            break;
    }

    configData.T = T;
    configData.Cm = Cm;
    configData.Ct = Ct;

    // Calculate diameter
    const shaftType = document.querySelector('input[name="shaftType"]:checked').value;
    const maxBM = configData.maxBM;
    const term1 = Cm * maxBM;
    const term2 = Math.sqrt(Math.pow(Cm * maxBM, 2) + Math.pow(Ct * T, 2));

    let diameter, outerDia, innerDia, dNormal, dShear, diaDetails = '';

    if (shaftType === 'solid') {
        dNormal = Math.pow((16 / (
        result = calcUDL(T);
    } else if (problemType === 'twogear') {
        result = calcTwoGear(T);
    } else if (problemType === 'pulleygear') {
        result = calcPulleyGear(T);
    }

    // Calculate shaft diameter
    const shaftType = document.querySelector('input[name="shaftType"]:checked').value;
    let diameter, outerDia, innerDia;
    const maxBM = result.maxBM;

    const term1 = Cm * maxBM;
    const term2 = Math.sqrt(Math.pow(Cm * maxBM, 2) + Math.pow(Ct * T, 2));

    const dNormal = Math.pow((16 / (Math.PI * sigmaMax)) * (term1 + term2), 1 / 3);
    const dShear = Math.pow((16 / (Math.PI * tauMax)) * term2, 1 / 3);

    if (shaftType === 'solid') {
        diameter = Math.max(dNormal, dShear);
        const stdSizes = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120];
        diameter = stdSizes.find(sMath.PI * stress.sigmaMax)) * (term1 + term2), 1 / 3);
        dShear = Math.pow((16 / (Math.PI * stress.tauMax)) * term2, 1 / 3);
        diameter = Math.max(dNormal, dShear);

        diaDetails += `<h3>Diameter Calculation</h3>`;
        diaDetails += `<p><strong>Method 1 - Max Normal Stress Theory:</strong></p>`;
        diaDetails += `<p>d = [16/(π×${stress.sigmaMax.toFixed(2)}) × {${Cm}×${maxBM.toFixed(0)} + √((${Cm}×${maxBM.toFixed(0)})² + (${Ct}×${T.toFixed(0)})²)}]^(1/3)</p>`;
        diaDetails += `<p><strong>d = ${dNormal.toFixed(2)} mm ... Eq.(i)</strong></p>`;
        diaDetails += `<p><strong>Method 2 - Max Shear Stress Theory:</strong></p>`;
        diaDetails += `<p>d = [16/(π×${stress.tauMax.toFixed(2)}) × √((${Cm}×${maxBM.toFixed(0)})² + (${Ct}×${T.toFixed(0)})²)]^(1/3)</p>`;
        diaDetails += `<p><strong>d = ${dShear.toFixed(2)} mm ... Eq.(ii)</strong></p> => s >= diameter) || Math.ceil(diameter / 10) * 10;
        outerDia = diameter;
        innerDia = 0;
    } else {
        const k = evalFraction(document.getElementById('diameterRatio').value);
        const doNormal = Math.pow((16 / (Math.PI * sigmaMax)) * (term1 + term2) * (1 / (1 - Math.pow(k, 4))), 1 / 3);
        const doShear = Math.pow((16 / (Math.PI * tauMax)) * term2 * (1 / (1 - Math.pow(k, 4))), 1 / 3);
        outerDia = Math.max(doNormal, doShear);
        const stdSizes = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120];
        outerDia = stdSizes.find(s => s >= outerDia) || Math.ceil(outerDia / 10) * 10;
        innerDia = k * outerDia;
        diameter = outerDia;
    }

    // Display results
    displayResults(T, maxBM, diameter, outerDia, innerDia, shaftType, Cm, Ct, sigmaMax, tauMax, dNormal, dShear, result);
    generateDiagrams(result, T);
}

// ========================================
// POINT LOAD CALCULATION
// ========`;
        diaDetails += `<p><strong>Select maximum: d = ${diameter.toFixed(2)} mm</strong></p>`;

        const stds = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120];
        outerDia = stds.find(s => s >= diameter) || Math.ceil(diameter / 10) * 10;
        innerDia = 0;
        diaDetails += `<p><strong>∴ Standard size: d = ${outerDia} mm</strong></p>`;
    } else {
        const kStr = document.getElementById('diameterRatio').value;
        const k = evalFraction(kStr);
        const kFactor = 1 / (1 - Math.pow(k, 4));

        dNormal = Math.pow((16 / (Math.PI * stress.sigmaMax)) * (term1 + term2) * kFactor, 1 / 3);
        dShear = Math.pow((16 / (Math.PI * stress.tauMax)) * term2 * kFactor, 1 / 3);
        outerDia = Math.max(dNormal, dShear);

        const stds = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120];
        outerDia = stds.find================================
function calcPointLoad(T) {
    const L = 600;
    const b1 = 100, b2 = 500, gp = 300;
    const Dg = 200;
    const Ft = (2 * T) / Dg;
    const Fr = Ft * Math.tan(20 * Math.PI / 180);

    const RAh = Ft * (b2 - gp) / (b2 - b1);
    const RBh = Ft - RAh;
    const RAv = Fr * (b2 - gp) / (b2 - b1);
    const RBv = Fr - RAv;

    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        if (x[i] < b1) { SF_h.push(0); BM_h.push(0); SF_v.push(0); BM_v.push(0); }
        else if (x[i] < gp) {
            SF_h.push(RAh); BM_h.push(RAh * (x[i] - b1));
            SF_v.push((s => s >= outerDia) || Math.ceil(outerDia / 10) * 10;
        innerDia = k * outerDia;
        diameter = outerDia;

        diaDetails += `<h3>Diameter Calculation (Hollow Shaft, k=${k.toFixed(4)})</h3>`;
        diaDetails += `<p>Normal Stress: do = ${dNormal.toFixed(2)} mm ... Eq.(i)</p>`;
        diaDetails += `<p>Shear Stress: do = ${dShear.toFixed(2)} mm ... Eq.(ii)</p>`;
        diaDetails += `<p><strong>∴ Standard: do = ${outerDia} mm, di = ${innerDia.toFixed(2)} mm</strong></p>`;
    }

    // Display results
    document.getElementById('results').style.display = 'block';
    document.getElementById('diagrams').style.display = 'block';

    document.getElementById('torqueResult').textContent = `${T.toFixed(2)} N-mm`;
    document.getElementById('momentResult').textContent = `${maxBM.toFixed(2)} N-mm`;

    if (shaftType === 'solid') {
        document.getElementById('diameterResult').textContent = `d = ${outerDia} mm`;
    } else {
        document.getElementById('diameterResult').textContent = `do=${outerDia} mm, di=${innerDia.toFixed(2)} mm`;
    }

    let detHTML = `<h3>Step-by-StepRAv); BM_v.push(RAv * (x[i] - b1));
        } else if (x[i] <= b2) {
            SF_h.push(RAh - Ft); BM_h.push(RAh * (x[i] - b1) - Ft * (x[i] - gp));
            SF_v.push(RAv - Fr); BM_v.push(RAv * (x[i] - b1) - Fr * (x[i] - gp));
        } else { SF_h.push(0); BM_h.push(0); SF_v.push(0); BM_v.push(0); }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return { x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM, type: 'point', bearing1: b1, bearing2: b2 };
}

// ========================================
// UDL CALCULATION
// ========================================
function calcUDL(T) {
    const L = parseFloat(document.getElementById('totalLength').value);
    const lUdl = parseFloat(document.getElementById('udlLength').value);
    const w = parseFloat(document.getElementById('udlIntensity').value);

    const b1 = (L - lUdl) / 2;
    const b2 = b1 + lUdl;
    const RA = (w * lUdl) / 2;
    const n = 500;
    const x Solution</h3>`;
    detHTML += `<p><strong>Torque:</strong> T = 9.55×10⁶×${P}/${N} = ${T.toFixed(2)} N-mm</p>`;
    detHTML += configData.details || '';
    detHTML += `<h3>Material Stress Analysis</h3>`;
    detHTML += stress.details;
    detHTML += `<p>Cm = ${Cm}, Ct = ${Ct}</p>`;
    detHTML += diaDetails;
    document.getElementById('detailedResults').innerHTML = detHTML;

    // Generate diagrams
    generateDiagrams(configData);

    document.getElementById('results').scrollIntoView({ behavior: 'smooth' });
}

// ========== POINT LOAD ==========
function calcPointLoad(T) {
    const L = 600, b1 = 100, b2 = 500, gp = 300;
    const Dg = 200;
    const Ft = (2 * T) / Dg;
    const Fr = Ft * Math.tan(20 * Math.PI / 180);

    const RAh = Ft * (b2 - gp) / (b2 - b1);
    const RBh = Ft - RAh;
    const RAv = Fr * = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        if (x[i] < b1 || x[i] > b2) {
            SF_h.push(0); BM_h.push(0); SF_v.push(0); BM_v.push(0);
        } else {
            const d = x[i] - b1;
            const sf = RA - w * d;
            const bm = RA * d - w * d * d / 2;
            SF_h.push(sf); BM_h.push(bm); SF_v.push(sf); BM_v.push(bm);
        }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return { x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM, type: 'udl', bearing1: b1, bearing2: b2 };
}

// ========================================
// TWO GEAR CALCULATION (Problems 21, 22)
// ========================================
function calcTwoGear(T) {
    const AB = parseFloat(document.getElementById('bearingDistance').value);
    const AC = parseFloat(document.getElementById('gearC_pos').value);
    const AD = parseFloat(document.getElementById('gearD_pos').value);
    const Dc = parseFloat(document.getElementById('gearC_dia').value);
    const Dd = parseFloat(document.getElementById('gearD_dia').value);
    const W (b2 - gp) / (b2 - b1);
    const RBv = Fr - RAv;

    const BMh = RAh * (gp - b1);
    const BMv = RAv * (gp - b1);
    const maxBM = Math.sqrt(BMh * BMh + BMv * BMv);

    return { type: 'point', L, b1, b2, gp, Ft, Fr, RAh, RBh, RAv, RBv, maxBM };
}

// ========== UDL ==========
function calcUDL(T) {
    const L = parseFloat(document.getElementById('totalLength').value);
    const lUdl = parseFloat(document.getElementById('udlLength').value);
    const w = parseFloat(document.getElementById('udlIntensity').value);
    const b1 = (L - lUdl) / 2;
    const b2 = b1 + lUdl;
    const WTotal = w * lUdl;
    const RA = WTotal / 2;
    const RB = WTotal / 2;
    const maxBM = (w * lUdl * lUdl) / 8;

    let details = `<h3>UDL Analysis</h3>`;
    details += `<p>Total UDL force: W = ${w} × ${lUdl} = ${WTotal} N</p>`;
    details += `<p>RA = RB = ${RA.toFixed(2)} N</p>`;
    details += `<p>Max Bc = parseFloat(document.getElementById('gearC_weight').value);
    const Wd = parseFloat(document.getElementById('gearD_weight').value);
    const alpha = parseFloat(document.getElementById('twoGear_pressureAngle').value) * Math.PI / 180;

    const FtC = (2 * T) / Dc;
    const FrC = FtC * Math.tan(alpha);
    const FtD = (2 * T) / Dd;
    const FrD = FtD * Math.tan(alpha);

    const L = AB + 50;

    // VERTICAL PLANE: Forces are FrD (down at D), FrC (up at C), weights
    // Taking moments about A
    const RBv = (FtD * AD - FrC * AC + Wc * AC + (FtD > FrC ? Wd * AD : Wd * AD)) / AB;
    // Simplified: use standard moment equation
    const RBv_calc = ((FrD + Wd) * AD - (FrC - Wc) * AC) / AB;
    const RAv_calc = (FrD + Wd) + (FrC +M = wl²/8 = ${maxBM.toFixed(2)} N-mm</p>`;

    return { type: 'udl', L, b1, b2, w, lUdl, RA, RB, maxBM, details };
}

// ========== TWO GEAR ==========
function calcTwoGear(T) {
    const AB = parseFloat(document.getElementById('bearingDist').value);
    const AC = parseFloat(document.getElementById('gearCPos').value);
    const AD = parseFloat(document.getElementById('gearDPos').value);
    const Dc = parseFloat(document.getElementById('gearCDia').value);
    const Dd = parseFloat(document.getElementById('gearDDia').value);
    const alpha = parseFloat(document.getElementById('pressureAngle2g').value) * Math.PI / 180;
    const Wc = parseFloat(document.getElementById('gearCWeight').value);
    const Wd = parseFloat(document.getElementById('gearDWeight').value);

    const FtC = (2 * T) / Dc;
    const FrC = FtC * Math.tan(alpha);
    const FtD = (2 * T) / Dd;
    const FrD = FtD * Math.tan(alpha);

    // Vertical plane: radial forces + weights act vertically
    const FvC = FrC + Wc;
    const FvD = FrD + Wd;

    // Vertical: Taking moments about A
    const RBV Wc) - RBv_calc;

    // Recalculate properly based on textbook method
    // Vertical: FtD acts at D, FrC acts at C (both can be up or down)
    // For simplicity, use absolute values and directions from textbook
    const RBv2 = (FtD * AD - FrC * AC) / AB;
    const RAv2 = FtD - FrC - RBv2;

    // HORIZONTAL PLANE
    const RBh = (FtC * AC - FtD * AD) / AB;
    const RAh_abs = Math.abs(RBh) + FtC - FtD;

    // Use proper method from textbook
    // Vertical plane reactions
    const V_RB = (FrD * AD - FrC * AC) / AB;
    const V_RA = FrD - FrC + V_RB;

    // Horizontal plane reactions
    const H_RB = (FtC * AC + FtD * AD) / AB;
    const H_RA = FtC + FtD - H_RB;

    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        // Vertical plane
        if (x = (FtD * AD - FrC * AC) / AB;
    const RAV = FtD - FrC - RBV;

    // Actually use proper method from textbook
    // Vertical: FrC at C, FrD (and FtD) at D
    const RBv = (FrD * AD - FrC * AC) / AB;
    const RAv = FrD - FrC - RBv;

    // Horizontal: FtC at C, FtD at D
    const RBh_raw = (FtC * AC - FtD * AD) / AB;
    const RAh_raw = FtC - FtD - RBh_raw;

    // Taking moments about A for horizontal
    // RBH × AB + FtD × AD = FtC × AC
    const RBh = (FtC * AC - FtD * AD) / AB;
    const RAh = FtC + FtD + RBh; // Force balance

    const L = AB;
    const b1 = 0;
    const b2 = AB;

    // Vertical BM at key points
    const MCV = RAv * AC;
    const MDV = RBv * (AB[i] < 0) { SF_v.push(0); BM_v.push(0); }
        else if (x[i] < AC) { SF_v.push(V_RA); BM_v.push(V_RA * x[i]); }
        else if (x[i] < AD) { SF_v.push(V_RA - FrC); BM_v.push(V_RA * x[i] - FrC * (x[i] - AC)); }
        else if (x[i] <= AB) { SF_v.push(V_RA - FrC + FrD); BM_v.push(V_RA * x[i] - FrC * (x[i] - AC) + FrD * (x[i] - AD)); }
        else { SF_v.push(0); BM_v.push(0); }

        // Horizontal plane
        if (x[i] < 0) { SF_h.push(0); BM_h.push(0); }
        else if (x[i] < AC) { SF_h.push(H_RA); BM_h.push(H_RA * x[i]); }
        else if (x[i] < AD) { SF_h.push(H_RA - FtC); BM_h.push(H_RA * x[i] - FtC * (x[i] - AC)); }
        else if (x[i] <= AB) { SF_h.push(H_RA - FtC + FtD); BM_h.push(H_RA * x[i] - FtC * (x[i] - AC) + FtD * (x[i] - AD)); }
        else { SF_h.push(0); BM_h.push(0); }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'twog - AD);

    // Horizontal BM at key points
    const MCH = RAh * AC;
    const MDH = -RBh * (AB - AD);

    // Resultant moments
    const MC = Math.sqrt(MCV * MCV + MCH * MCH);
    const MD = Math.sqrt(MDV * MDV + MDH * MDH);
    const maxBM = Math.max(MC, MD);

    let details = `<h3>Two Gear Analysis</h3>`;
    details += `<p><strong>Gear C:</strong> Ft_C = 2T/D_C = ${FtC.toFixed(2)} N, Fr_C = Ft_C×tan(α) = ${FrC.toFixed(2)} N</p>`;
    details += `<p><strong>Gear D:</strong> Ft_D = 2T/D_D = ${FtD.toFixed(2)} N, Fr_D = Ft_D×tan(α) = ${FrD.toFixed(2)} N</p>`;
    details += `<h3>Vertical Loading</h3>`;
    details += `<p>R_AV = ${Math.abs(RAv).toFixed(2)} N, R_BV = ${Math.abs(RBv).toFixed(2)} N</p>`;
    details += `<p>M_CV = ${Math.abs(MCV).toFixed(2)} N-mm, M_DV = ${Math.abs(MDV).toFixed(2)} N-mm</p>`;
    details += `<h3>Horizontal Loading</h3>`;
    details += `<p>R_AH = ${Math.abs(RAh).toFixed(2)} N, R_BH = ${Math.abs(RBh).toFixed(2)} N</p>`;
    details += `<p>M_CH = ${Math.abs(MCear', bearing1: 0, bearing2: AB,
        FtC, FrC, FtD, FrD, AC, AD, AB,
        V_RA, V_RB, H_RA, H_RB
    };
}

// ========================================
// PULLEY + GEAR CALCULATION (Problems 23, 24)
// ========================================
function calcPulleyGear(T) {
    const CD = parseFloat(document.getElementById('pg_bearingDistance').value);
    const pulleyPos = parseFloat(document.getElementById('pulley_pos').value);
    const pulleyDia = parseFloat(document.getElementById('pulley_dia').value);
    const Wa = parseFloat(document.getElementById('pulley_weight').value);
    const tensionRatio = parseFloat(document.getElementById('tensionRatio').value);
    const gearPos = parseFloat(document.getElementById('gear_pos').value);
    const gearDia = parseFloat(document.getElementById('gear_dia').value);
    const Wb = parseFloat(document.getElementById('gear_weight').value);
    const alpha = parseFloat(document.getElementById('pg_pressureAngle').value) * Math.PI / 180;

    // Pulley: T_A = (T1 - T2) * R_pulley
    // T1/T2 = tensionRatio
    const R_pulley = pulleyDia / 2;
    const T2 = T / (R_pulley * (tensH).toFixed(2)} N-mm, M_DH = ${Math.abs(MDH).toFixed(2)} N-mm</p>`;
    details += `<h3>Resultant Bending Moments</h3>`;
    details += `<p>M_C = √(M_CV² + M_CH²) = ${MC.toFixed(2)} N-mm</p>`;
    details += `<p>M_D = √(M_DV² + M_DH²) = ${MD.toFixed(2)} N-mm</p>`;
    details += `<p><strong>Max BM = ${maxBM.toFixed(2)} N-mm</strong></p>`;

    return {
        type: 'twogear', L, b1, b2, AC, AD,
        FtC, FrC, FtD, FrD,
        RAv, RBv, RAh, RBh,
        MCV, MDV, MCH, MDH, MC, MD,
        maxBM, details
    };
}

// ========== PULLEY + GEAR ==========
function calcPulleyGear(T) {
    const CD = parseFloat(document.getElementById('pgBearingDist').value);
    const pulleyFromC = parseFloat(document.getElementById('pulleyPos').value);
    const gearFromC = parseFloat(document.getElementById('gearPosAdv').value);
    const Dp = parseFloat(document.getElementById('pulleyDia').value);
    const Dg = parseFloat(document.getElementById('gearDiaAdv').value);
    const T1 = parseFloat(document.getElementById('tension1').value);
    const T2 = parseFloat(document.getElementById('tension2').value);
    const Wp = parseFloat(document.getElementById('pulleyWeight').value);
    const Wg = parseFloat(ionRatio - 1));
    const T1 = tensionRatio * T2;

    // Gear forces
    const FtB = (2 * T) / gearDia;
    const FrB = FtB * Math.tan(alpha);

    // Layout: A---C---D---B (C=pulley, D=gear, bearings at C and D)
    const AC = pulleyPos;
    const DB = gearPos - CD;
    const L = pulleyPos + CD + DB + 50;

    // Vertical plane
    // At pulley: T1 + T2 + Wa acts downward
    // At gear: FtB - Wb (or FtB + Wb depending on direction)
    const vertPulleyForce = T1 + T2 + Wa;
    const vertGearForce = FtB + Wb;

    // Reactions at bearings C and D (vertical)
    const RDv = (vertPulleyForce * AC + vertGearForce * (CD + DB)) / CD;
    const RCv = vertPulleyForce + vertGearForce - RDv;

    // Horizontal plane
    const horizGearForce = FrB;
    const RDh = horizGearForce * DB / CD;
    const RCh = horizGearForce -document.getElementById('gearWeightAdv').value);
    const alpha = parseFloat(document.getElementById('pressureAnglePG').value) * Math.PI / 180;

    // Gear forces
    const FtB = (2 * T) / Dg;
    const FrB = FtB * Math.tan(alpha);

    // Pulley: belt tensions act vertically
    const Fv_pulley = T1 + T2 + Wp;
    const Fv_gear = FrB + Wg;

    // For vertical loading
    const AC = Math.abs(pulleyFromC);
    const CB = gearFromC > CD ? CD : gearFromC;

    // Vertical: moments about C
    let RDV, RCV;
    if (pulleyFromC < 0) {
        // Pulley is to the left of C
        RDV = ((FtB + Wg) * gearFromC - (T1 + T2 + Wp) * Math.abs(pulleyFromC)) / CD;
        RCV = (T1 + T2 + Wp) + (FtB + Wg) - RDV;
    } else {
        RDV = ((T1 + T2 + Wp) * pulleyFromC + (FtB + Wg) * gearFromC) / CD;
        RCV = (T1 + T2 + Wp) + (FtB + Wg) - RDV;
    }

    // Horizontal: moments about C
    let RDH, RCH; RDh;

    const n = 500;
    const bearingC = AC;
    const bearingD = AC + CD;
    const gearPosAbs = AC + CD + DB;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        // Vertical plane
        if (x[i] < 0) { SF_v.push(0); BM_v.push(0); }
        else if (x[i] < bearingC) {
            SF_v.push(-vertPulleyForce);
            BM_v.push(-vertPulleyForce * x[i]);
        } else if (x[i] < bearingD) {
            SF_v.push(-vertPulleyForce + RCv);
            BM_v.push(-vertPulleyForce * x[i] + RCv * (x[i] - bearingC));
        } else if (x[i] <= gearPosAbs) {
            SF_v.push(-vertPulleyForce + RCv + RDv);
            BM_v.push(-vertPulleyForce * x[i] + RCv * (x[i] - bearingC) + RDv * (x[i] - bearingD));
        } else { SF_v.push(0); BM_v.push(0); }

        // Horizontal plane
        if (x[i] < bearingC) { SF_h.push(0); BM_h.push(0); }
        else if (x[i] < bearingD) {
            SF_h.push(RCh);
            BM_h.push(RCh * (x[i] - bearingC));
        } else if (x[i] <= gearPosAbs) {
            SF_h.push(RCh - RDh);
            BM
    if (gearFromC > CD) {
        RDH = (FrB * gearFromC) / CD;
        RCH = FrB - RDH;
    } else {
        RDH = (FrB * gearFromC) / CD;
        RCH = FrB - RDH;
    }

    // Bending moments at key points
    let MCV, MDV, MCH, MDH;

    if (pulleyFromC < 0) {
        MCV = (T1 + T2 + Wp) * Math.abs(pulleyFromC);
    } else {
        MCV = RCV * pulleyFromC;
    }

    MDV = RDV * (CD - gearFromC > 0 ? CD - gearFromC : gearFromC - CD);
    if (gearFromC > CD) {
        MDV = (FtB + Wg) * (gearFromC - CD);
    }

    MCH = 0;
    if (gearFromC > CD) {
        MDH = FrB * (gearFromC - CD);
    } else {
        MDH = RDH * (CD - gearFromC);
    }

    const MC = Math.sqrt(MCV * MCV + MCH * MCH);
    const MD = Math.sqrt(MDV * MDV + MDH * MDH);
    const maxBM = Math.max(MC, MD);

    const L = Math.max(CD, gearFromC, Math.abs(pulleyFromC)) + 100;

    let details = `<h3>Pulley + Gear Analysis</h3>`;
    details += `<p><strong>Gear:</strong> Ft =_h.push(RCh * (x[i] - bearingC) + RDh * (x[i] - bearingD));
        } else { SF_h.push(0); BM_h.push(0); }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'pulleygear', bearing1: bearingC, bearing2: bearingD,
        T1, T2, FtB, FrB, Wa, Wb, RCv, RDv, RCh, RDh,
        vertPulleyForce, vertGearForce, horizGearForce
    };
}

// ========================================
// DISPLAY RESULTS
// ========================================
function displayResults(T, maxBM, diameter, outerDia, innerDia, shaftType, Cm, Ct, sigmaMax, tauMax, dNormal, dShear, result) {
    document.getElementById('results').style.display = 'block';
    document.getElementById('diagrams').style.display = 'block';

    document.getElementById('torqueResult').textContent = T.toFixed(2) + ' N-mm (' + (T / 1000).toFixed(2) + ' N-m)';
    document.getElementById('momentResult').textContent = maxBM.toFixed(2) + ' N-mm';

    if (shaftType === 'solid') {
        document.getElementById('diameterResult').textContent = 'd = ' + diameter.toFixed(2 ${FtB.toFixed(2)} N, Fr = ${FrB.toFixed(2)} N</p>`;
    details += `<p><strong>Pulley:</strong> T1 = ${T1} N, T2 = ${T2} N</p>`;
    details += `<h3>Vertical Loading</h3>`;
    details += `<p>R_CV = ${Math.abs(RCV).toFixed(2)} N, R_DV = ${Math.abs(RDV).toFixed(2)} N</p>`;
    details += `<p>M_CV = ${Math.abs(MCV).toFixed(2)} N-mm</p>`;
    details += `<p>M_DV = ${Math.abs(MDV).toFixed(2)} N-mm</p>`;
    details += `<h3>Horizontal Loading</h3>`;
    details += `<p>R_CH = ${Math.abs(RCH).toFixed(2)} N, R_DH = ${Math.abs(RDH).toFixed(2)} N</p>`;
    details += `<p>M_DH = ${Math.abs(MDH).toFixed(2)} N-mm</p>`;
    details += `<h3>Resultant Bending Moments</h3>`;
    details += `<p>M_C = ${MC.toFixed(2)} N-mm</p>`;
    details += `<p>M_D = ${MD.toFixed(2)} N-mm</p>`;
    details += `<p><strong>Max BM = ${maxBM.toFixed(2)} N-mm</strong></p>`;

    return {
        type: 'pulleygear', L, b1: 0, b2: CD, CD,
        pulleyFromC, gearFromC,
        FtB, FrB, T1, T2, Wp, Wg,
        RCV, RDV, RCH, RDH,
        MCV, MDV, MCH, MDH, MC, MD,) + ' mm';
    } else {
        document.getElementById('diameterResult').textContent = 'do = ' + outerDia.toFixed(2) + ' mm, di = ' + innerDia.toFixed(2) + ' mm';
    }

    let detailsHTML = '<h3>Step-by-Step Calculation</h3>';
    detailsHTML += '<p><strong>Torque:</strong> T = 9.55×10⁶ × P / n = ' + T.toFixed(2) + ' N-mm</p>';
    detailsHTML += '<p><strong>σ_max:</strong> ' + sigmaMax.toFixed(2) + ' MPa</p>';
    detailsHTML += '<p><strong>τ_max:</strong> ' + tauMax.toFixed(2) + ' MPa</p>';
    detailsHTML += '<p><strong>Cm = ' + Cm.toFixed(2) + ', Ct = ' + Ct.toFixed(2) + '</strong></p>';
    detailsHTML += '<p><strong>Maximum Bending Moment M:</strong> ' + maxBM.toFixed(2) + ' N-mm</p>';

    if (result.type === 'twogear') {
        detailsHTML += '<h3>Gear Force Analysis</h3>';
        detailsHTML += '<p>Gear C: Ft = ' + result.FtC.toFixed(2) + ' N, Fr = ' + result.FrC.toFixed(2) + ' N</p>';
        detail
        maxBM, details
    };
}

// ========== GENERATE DIAGRAMS ==========
function generateDiagrams(cfg) {
    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * cfg.L);

    let SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    if (cfg.type === 'point') {
        for (let i = 0; i < n; i++) {
            if (x[i] < cfg.b1) { SF_h.push(0); BM_h.push(0); SF_v.push(0); BM_v.push(0); }
            else if (x[i] < cfg.gp) {
                SF_h.push(cfg.RAh); BM_h.push(cfg.RAh * (x[i] - cfg.b1));
                SF_v.push(cfg.RAv); BM_v.push(cfg.RAv * (x[i] - cfg.b1));
            } else if (x[i] <= cfg.b2) {
                SF_h.push(cfg.RAh - cfg.Ft); BM_h.push(cfg.RAh * (x[i] - cfg.b1) - cfg.Ft * (x[i] - cfg.gp));
                SF_v.push(cfg.RAv - cfg.Fr); BM_v.push(cfg.RAv * (x[i] - cfg.b1) - cfg.Fr * (x[i] - cfg.gp));
            } else { SF_h.push(0); BM_h.push(0); SF_v.push(0); BM_v.push(0); }
        }sHTML += '<p>Gear D: Ft = ' + result.FtD.toFixed(2) + ' N, Fr = ' + result.FrD.toFixed(2) + ' N</p>';
        detailsHTML += '<h3>Reactions</h3>';
        detailsHTML += '<p>Vertical: RA = ' + result.V_RA.toFixed(2) + ' N, RB = ' + result.V_RB.toFixed(2) + ' N</p>';
        detailsHTML += '<p>Horizontal: RA = ' + result.H_RA.toFixed(2) + ' N, RB = ' + result.H_RB.toFixed(2) + ' N</p>';
    }

    if (result.type === 'pulleygear') {
        detailsHTML += '<h3>Force Analysis</h3>';
        detailsHTML += '<p>Pulley: T1 = ' + result.T1.toFixed(2) + ' N, T2 = ' + result.T2.toFixed(2) + ' N</p>';
        detailsHTML += '<p>Gear: Ft = ' + result.FtB.toFixed(2) + ' N, Fr = ' + result.FrB.toFixed(2) + ' N</p>';
        detailsHTML += '<h3>Reactions</h3>';
        detailsHTML += '<p>Vertical: RC = ' + result.RCv.toFixed(2) + ' N, RD = ' + result.RDv.toFixed(2) + ' N</p>';
        detailsHTML += '<p>Horizontal: RC = ' + result.RCh.toFixed(2) + ' N, RD = ' + result.RDh.toFixed(2) + ' N</p>';
    }

    detailsHTML += '<h3>Diameter Calculation</h3>';
    detailsHTML += '<p><strong>Eq.(i) Normal Stress Theory
    } else if (cfg.type === 'udl') {
        for (let i = 0; i < n; i++) {
            if (x[i] < cfg.b1 || x[i] > cfg.b2) { SF_h.push(0); BM_h.push(0); SF_v.push(0); BM_v.push(0); }
            else {
                const d = x[i] - cfg.b1;
                const sf = cfg.RA - cfg.w * d;
                const bm = cfg.RA * d - cfg.w * d * d / 2;
                SF_h.push(sf); BM_h.push(bm); SF_v.push(sf); BM_v.push(bm);
            }
        }
    } else if (cfg.type === 'twogear') {
        for (let i = 0; i < n; i++) {
            // Vertical
            if (x[i] < cfg.b1) { SF_v.push(0); BM_v.push(0); }
            else if (x[i] < cfg.AC) { SF_v.push(cfg.RAv); BM_v.push(cfg.RAv * x[i]); }
            else if (x[i] < cfg.AD) { SF_v.push(cfg.RAv - cfg.FrC); BM_v.push(cfg.RAv * x[i] - cfg.FrC * (x[i] - cfg.AC)); }
            else if (x[i] <= cfg.b2) { SF_v.push(cfg.RAv - cfg.FrC + cfg.FrD); BM_v.push(cfg.RAv * x[i] - cfg.FrC * (x[i] - cfg.AC) + cfg.FrD * (x[i] - cfg.AD)); }
            else { SF_v.push(0); BM_v.push(0); }:</strong> d = ' + dNormal.toFixed(2) + ' mm</p>';
    detailsHTML += '<p><strong>Eq.(ii) Shear Stress Theory:</strong> d = ' + dShear.toFixed(2) + ' mm</p>';
    detailsHTML += '<p><strong>Design Diameter:</strong> d = ' + Math.max(dNormal, dShear).toFixed(2) + ' mm</p>';

    if (shaftType === 'solid') {
        detailsHTML += '<p style="color:red; font-weight:bold;">∴ Standard size of shaft: d = ' + outerDia + ' mm</p>';
    } else {
        detailsHTML += '<p style="color:red; font-weight:bold;">∴ Standard size: do = ' + outerDia + ' mm, di = ' + innerDia.toFixed(2) + ' mm</p>';
    }

    document.getElementById('detailedResults').innerHTML = detailsHTML;
    document.getElementById('results').scrollIntoView({ behavior: 'smooth' });
}

// ========================================
// GENERATE DIAGRAMS
// ========================================
function generateDiagrams(result, T) {
    const x = result.x;
    const n = x.length;

    const torqueData = result.type === 'udl' ? Array(n).fill(T) :
        (result.type === 'point' ? x.map(xi => xi >= 300 ? T : 0) :

            // Horizontal
            if (x[i] < cfg.b1) { SF_h.push(0); BM_h.push(0); }
            else if (x[i] < cfg.AC) { SF_h.push(cfg.RAh); BM_h.push(cfg.RAh * x[i]); }
            else if (x[i] < cfg.AD) { SF_h.push(cfg.RAh - cfg.FtC); BM_h.push(cfg.RAh * x[i] - cfg.FtC * (x[i] - cfg.AC)); }
            else if (x[i] <= cfg.b2) { SF_h.push(cfg.RAh - cfg.FtC + cfg.FtD); BM_h.push(cfg.RAh * x[i] - cfg.FtC * (x[i] - cfg.AC) + cfg.FtD * (x[i] - cfg.AD)); }
            else { SF_h.push(0); BM_h.push(0); }
        }
    } else if (cfg.type === 'pulleygear') {
        const pPos = cfg.pulleyFromC < 0 ? 0 : cfg.pulleyFromC;
        const gPos = cfg.gearFromC > cfg.CD ? cfg.CD : cfg.gearFromC;

        for (let i = 0; i < n; i++) {
            // Simplified vertical
            if (x[i] < 0) { SF_v.push(0); BM_v.push(0); }
            else if (x[i] <= cfg.CD) { SF_v.push(cfg.RCV - (x[i] > pPos ? cfg.T1 + cfg.T2 + cfg.Wp : 0)); BM_v.push(cfg.RCV * x[i] - (x[i] > pPos ? (cfg.T1 + cfg.T2 + cfg.Wp) * (x[i] - pPos) : 0)); } Array(n).fill(T));

    document.getElementById('diagrams').innerHTML = '<h2>📈 Horizontal Loading Diagram (HLD)</h2>' +
        '<div id="sfd_h"></div><div id="bmd_h"></div>' +
        '<h2>📈 Vertical Loading Diagram (VLD)</h2>' +
        '<div id="sfd_v"></div><div id="bmd_v"></div>' +
        '<h2>📈 Resultant Bending Moment</h2><div id="bmd_res"></div>' +
        '<h2>📈 Torque Diagram</h2><div id="torque_d"></div>';

    // HLD - SFD
    Plotly.newPlot('sfd_h', [{ x: x, y: result.SF_h, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(55,128,191)', width: 3 } }],
        { title: 'Shear Force Diagram - Horizontal Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'SF (N)' } });

    // HLD - BMD
    const maxH = Math.max(...result.BM_h.map(Math.abs));
    const idxH = result.BM_h.findIndex(b
            else { SF_v.push(0); BM_v.push(0); }

            // Simplified horizontal
            if (x[i] < 0) { SF_h.push(0); BM_h.push(0); }
            else if (x[i] <= cfg.CD) { SF_h.push(cfg.RCH - (x[i] > gPos ? cfg.FrB : 0)); BM_h.push(cfg.RCH * x[i] - (x[i] > gPos ? cfg.FrB * (x[i] - gPos) : 0)); }
            else { SF_h.push(0); BM_h.push(0); }
        }
    }

    // Resultant BM
    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const torqueData = cfg.type === 'udl' ? Array(n).fill(cfg.T) : (cfg.type === 'point' ? x.map(xi => xi >= cfg.gp ? cfg.T : 0) : Array(n).fill(cfg.T));

    // Build diagram HTML
    document.getElementById('diagrams').innerHTML = `
        <h2>📈 Horizontal Loading Diagram (HLD)</h2>
        <div id="sfd_h"></div><div id="bmd_h"></div>
        <h2>📈 Vertical Loading Diagram (VLD)</h => Math.abs(b) === maxH);
    Plotly.newPlot('bmd_h', [{ x: x, y: result.BM_h, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(219,64,82)', width: 3 } }],
        { title: 'Bending Moment Diagram - Horizontal Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'BM (N-mm)' },
            annotations: [{ x: x[idxH], y: result.BM_h[idxH], text: 'M_max=' + result.BM_h[idxH].toFixed(0), showarrow: true, arrowhead: 2, ax: 0, ay: -40 }] });

    // VLD - SFD
    Plotly.newPlot('sfd_v', [{ x: x, y: result.SF_v, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(55,128,191)', width: 3 } }],
        { title: 'Shear Force Diagram - Vertical Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'SF (N)' } });

    // VLD - BMD
    const maxV = Math.max(...result.BM_v.map(Math.abs));
    const idxV = result.BM_v.findIndex(b => Math.abs(b) === maxV);
    Plotly.newPlot('bmd_v', [{ x: x, y: result.BM_v, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(219,64,82)', width: 3 } }],
        { title: 'Bending Moment Diagram - Vertical Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'BM (N-mm)' },
            annotations: [{ x: x[idxV], y: result.BM_v[idxV], text: 'M_max=' + result.B2>
        <div id="sfd_v"></div><div id="bmd_v"></div>
        <h2>📈 Resultant Bending Moment</h2>
        <div id="bmd_res"></div>
        <h2>📈 Torque Diagram</h2>
        <div id="torque_d"></div>
    `;

    const plotCfg = { responsive: true };

    Plotly.newPlot('sfd_h', [{ x, y: SF_h, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(55,128,191)', width: 3 } }],
        { title: 'SFD - Horizontal Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'SF (N)' } }, plotCfg);

    const mxH = Math.max(...BM_h.map(Math.abs));
    const mxHi = BM_h.findIndex(b => Math.abs(b) === mxH);
    Plotly.newPlot('bmd_h', [{ x, y: BM_h, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(219,64,82)', width: 3 } }],
        { title: 'BMD - Horizontal Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'BM (N-mm)' }, annotations: [{ x: x[mxHi], y: BM_h[mxHi], text: `MM_v[idxV].toFixed(0), showarrow: true, arrowhead: 2, ax: 0, ay: -40 }] });

    // Resultant BMD
    const maxR = Math.max(...result.BM_res);
    const idxR = result.BM_res.indexOf(maxR);
    Plotly.newPlot('bmd_res', [{ x: x, y: result.BM_res, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(148,0,211)', width: 3 } }],
        { title: 'Resultant BM = √(M_h² + M_v²)', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'BM (N-mm)' },
            annotations: [{ x: x[idxR], y: maxR, text: 'M_max=' + maxR.toFixed(0), showarrow: true, arrowhead: 2, ax: 0, ay: -40 }] });

    // Torque
    Plotly.newPlot('torque_d', [{ x: x, y: torqueData, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(50,171,96)', width: 3 } }],
        { title: 'Torque Diagram', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'Torque (N-mm)' } });
}

document.addEventListener('DOMContentLoaded', function () { console.log('Shaft Design Calculator loaded'); });
