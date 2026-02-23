// ========================================
// TOGGLE FUNCTIONS
// ========================================
function toggleProblemType() {
    const val = document.querySelector('input[name="problemType"]:checked').value;
    document.getElementById('normalSection').style.display = val === 'normal' ? 'block' : 'none';
    document.getElementById('udlSection').style.display = val === 'udl' ? 'block' : 'none';
    document.getElementById('singlePulleySection').style.display = val === 'singlepulley' ? 'block' : 'none';
    document.getElementById('twoGearSection').style.display = val === 'twogear' ? 'block' : 'none';
    document.getElementById('pulleyGearSection').style.display = val === 'pulleygear' ? 'block' : 'none';
    document.getElementById('multiPulleySection').style.display = val === 'multipulley' ? 'block' : 'none';
}

function toggleStressMethod() {
    const method = document.querySelector('input[name="stressMethod"]:checked').value;
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

function toggleTensionInput() {
    const val = document.querySelector('input[name="tensionInput"]:checked').value;
    document.getElementById('directTensionInputs').style.display = val === 'direct' ? 'block' : 'none';
    document.getElementById('ratioTensionInputs').style.display = val === 'ratio' ? 'block' : 'none';
}

function toggleModulusG() {
    document.getElementById('modulusGInput').style.display = document.getElementById('hasModulusG').checked ? 'block' : 'none';
}

function toggleShaftLength() {
    document.getElementById('shaftLengthInput').style.display = document.getElementById('hasShaftLength').checked ? 'block' : 'none';
}

function toggleAngleTwist() {
    document.getElementById('angleTwistInput').style.display = document.getElementById('hasAngleTwist').checked ? 'block' : 'none';
}

function evalFraction(str) {
    str = str.trim();
    if (str.includes('/')) {
        const parts = str.split('/');
        return parseFloat(parts[0]) / parseFloat(parts[1]);
    }
    return parseFloat(str);
}

// ========================================
// MAIN CALCULATION
// ========================================
function calculateShaft() {
    const problemType = document.querySelector('input[name="problemType"]:checked').value;
    const P = parseFloat(document.getElementById('power').value);
    const N = parseFloat(document.getElementById('speed').value);
    const T = (9.55e6 * P) / N;

    // Get stress parameters
    const stressMethod = document.querySelector('input[name="stressMethod"]:checked').value;
    let sigmaMax, tauMax;

    if (stressMethod === 'direct') {
        tauMax = parseFloat(document.getElementById('tauAllow').value);
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

    // Keyway
    const hasKeyway = document.getElementById('hasKeyway').checked;
    if (hasKeyway) {
        sigmaMax = 0.75 * sigmaMax;
        tauMax = 0.75 * tauMax;
    }

    // Cm and Ct
    const useCmCt = document.getElementById('useCmCt').checked;
    const Cm = useCmCt ? parseFloat(document.getElementById('cm').value) : 1.5;
    const Ct = useCmCt ? parseFloat(document.getElementById('ct').value) : 1.0;

    // Calculate based on problem type
    let result;
    switch (problemType) {
        case 'normal': result = calcNormalShaft(T); break;
        case 'udl': result = calcUDL(T); break;
        case 'singlepulley': result = calcSinglePulley(T); break;
        case 'twogear': result = calcTwoGear(T); break;
        case 'pulleygear': result = calcPulleyGear(T); break;
        case 'multipulley': result = calcMultiPulley(T); break;
    }

    // Calculate diameter
    const shaftType = document.querySelector('input[name="shaftType"]:checked').value;
    const maxBM = result.maxBM;
    const term1 = Cm * maxBM;
    const term2 = Math.sqrt(Math.pow(Cm * maxBM, 2) + Math.pow(Ct * T, 2));

    const dNormal = Math.pow((16 / (Math.PI * sigmaMax)) * (term1 + term2), 1 / 3);
    const dShear = Math.pow((16 / (Math.PI * tauMax)) * term2, 1 / 3);

    let outerDia, innerDia, diameter;
    const stdSizes = [10, 12, 16, 20, 25, 30, 35, 40, 45, 50, 55, 60, 63, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120];

    if (shaftType === 'solid') {
        diameter = Math.max(dNormal, dShear);
        outerDia = stdSizes.find(s => s >= diameter) || Math.ceil(diameter / 10) * 10;
        innerDia = 0;
    } else {
        const k = evalFraction(document.getElementById('diameterRatio').value);
        const doNormal = Math.pow((16 / (Math.PI * sigmaMax)) * (term1 + term2) * (1 / (1 - Math.pow(k, 4))), 1 / 3);
        const doShear = Math.pow((16 / (Math.PI * tauMax)) * term2 * (1 / (1 - Math.pow(k, 4))), 1 / 3);
        outerDia = Math.max(doNormal, doShear);
        outerDia = stdSizes.find(s => s >= outerDia) || Math.ceil(outerDia / 10) * 10;
        innerDia = k * outerDia;
        diameter = outerDia;
    }

    displayResults(T, maxBM, diameter, outerDia, innerDia, shaftType, Cm, Ct, sigmaMax, tauMax, dNormal, dShear, result, hasKeyway);
    generateDiagrams(result, T);
}

// ========================================
// NORMAL SHAFT (Single Gear)
// ========================================
function calcNormalShaft(T) {
    const AB = parseFloat(document.getElementById('normal_bearingDist').value);
    const gearPos = parseFloat(document.getElementById('normal_gearPos').value);
    const gearDia = parseFloat(document.getElementById('normal_gearDia').value);
    const alpha = parseFloat(document.getElementById('normal_pressureAngle').value) * Math.PI / 180;
    const Wg = parseFloat(document.getElementById('normal_gearWeight').value);

    const Ft = (2 * T) / gearDia;
    const Fr = Ft * Math.tan(alpha);

    // Horizontal: Tangential force Ft at gear position
    const RAh = Ft * (AB - gearPos) / AB;
    const RBh = Ft - RAh;

    // Vertical: Radial force Fr + Weight at gear position
    const totalVertForce = Fr + Wg;
    const RAv = totalVertForce * (AB - gearPos) / AB;
    const RBv = totalVertForce - RAv;

    const L = AB;
    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        if (x[i] < gearPos) {
            SF_h.push(RAh); BM_h.push(RAh * x[i]);
            SF_v.push(RAv); BM_v.push(RAv * x[i]);
        } else {
            SF_h.push(RAh - Ft); BM_h.push(RAh * x[i] - Ft * (x[i] - gearPos));
            SF_v.push(RAv - totalVertForce); BM_v.push(RAv * x[i] - totalVertForce * (x[i] - gearPos));
        }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'normal', bearing1: 0, bearing2: AB,
        details: 'Ft = ' + Ft.toFixed(2) + ' N, Fr = ' + Fr.toFixed(2) + ' N\nRA_h = ' + RAh.toFixed(2) + ' N, RB_h = ' + RBh.toFixed(2) + ' N\nRA_v = ' + RAv.toFixed(2) + ' N, RB_v = ' + RBv.toFixed(2) + ' N',
        Ft, Fr, RAh, RBh, RAv, RBv
    };
}

// ========================================
// UDL SHAFT
// ========================================
function calcUDL(T) {
    const L = parseFloat(document.getElementById('totalLength').value);
    const lUdl = parseFloat(document.getElementById('udlLength').value);
    const w = parseFloat(document.getElementById('udlIntensity').value);

    const b1 = (L - lUdl) / 2;
    const b2 = b1 + lUdl;
    const RA = (w * lUdl) / 2;
    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
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

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'udl', bearing1: b1, bearing2: b2,
        details: 'w = ' + w + ' N/mm, l = ' + lUdl + ' mm\nRA = RB = ' + RA.toFixed(2) + ' N'
    };
}

// ========================================
// SINGLE PULLEY SHAFT
// ========================================
function calcSinglePulley(T) {
    const AB = parseFloat(document.getElementById('sp_bearingDist').value);
    const pulleyPos = parseFloat(document.getElementById('sp_pulleyPos').value);
    const pulleyDia = parseFloat(document.getElementById('sp_pulleyDia').value);
    const Wp = parseFloat(document.getElementById('sp_pulleyWeight').value);

    let T1, T2;
    const tensionInput = document.querySelector('input[name="tensionInput"]:checked').value;
    if (tensionInput === 'direct') {
        T1 = parseFloat(document.getElementById('sp_T1').value);
        T2 = parseFloat(document.getElementById('sp_T2').value);
    } else {
        const ratio = parseFloat(document.getElementById('sp_tensionRatio').value);
        const R_pulley = pulleyDia / 2;
        T2 = T / (R_pulley * (ratio - 1));
        T1 = ratio * T2;
    }

    // Bending moment M = (T1 + T2) * L at the pulley location
    const totalBeltForce = T1 + T2;
    const vertForce = totalBeltForce + Wp;

    // Reactions
    const RAv = vertForce * (AB - pulleyPos) / AB;
    const RBv = vertForce - RAv;

    // Horizontal: No horizontal force for vertical belt pull
    const RAh = 0;
    const RBh = 0;

    const L = AB;
    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        SF_h.push(0); BM_h.push(0);

        if (x[i] < pulleyPos) {
            SF_v.push(RAv); BM_v.push(RAv * x[i]);
        } else {
            SF_v.push(RAv - vertForce); BM_v.push(RAv * x[i] - vertForce * (x[i] - pulleyPos));
        }
    }

    const BM_res = BM_v.map(b => Math.abs(b));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'singlepulley', bearing1: 0, bearing2: AB,
        details: 'T1 = ' + T1.toFixed(2) + ' N, T2 = ' + T2.toFixed(2) + ' N\nBelt Force = ' + totalBeltForce.toFixed(2) + ' N\nRA = ' + RAv.toFixed(2) + ' N, RB = ' + RBv.toFixed(2) + ' N',
        T1, T2, RAv, RBv
    };
}

// ========================================
// TWO GEARS SHAFT
// ========================================
function calcTwoGear(T) {
    const AB = parseFloat(document.getElementById('bearingDistance').value);
    const AC = parseFloat(document.getElementById('gearC_pos').value);
    const AD = parseFloat(document.getElementById('gearD_pos').value);
    const Dc = parseFloat(document.getElementById('gearC_dia').value);
    const Dd = parseFloat(document.getElementById('gearD_dia').value);
    const Wc = parseFloat(document.getElementById('gearC_weight').value);
    const Wd = parseFloat(document.getElementById('gearD_weight').value);
    const alpha = parseFloat(document.getElementById('twoGear_pressureAngle').value) * Math.PI / 180;

    const FtC = (2 * T) / Dc;
    const FrC = FtC * Math.tan(alpha);
    const FtD = (2 * T) / Dd;
    const FrD = FtD * Math.tan(alpha);

    // Vertical plane: Radial forces (FrC, FrD) and weights
    const V_RB = ((FrD + Wd) * AD - (FrC + Wc) * AC) / AB;
    const V_RA = (FrD + Wd) - (FrC + Wc) + V_RB;

    // Horizontal plane: Tangential forces
    const H_RB = (FtC * AC + FtD * AD) / AB;
    const H_RA = FtC + FtD - H_RB;

    const L = AB;
    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        // Vertical
        if (x[i] < AC) {
            SF_v.push(V_RA); BM_v.push(V_RA * x[i]);
        } else if (x[i] < AD) {
            SF_v.push(V_RA - (FrC + Wc)); BM_v.push(V_RA * x[i] - (FrC + Wc) * (x[i] - AC));
        } else {
            SF_v.push(V_RA - (FrC + Wc) + (FrD + Wd)); BM_v.push(V_RA * x[i] - (FrC + Wc) * (x[i] - AC) + (FrD + Wd) * (x[i] - AD));
        }

        // Horizontal
        if (x[i] < AC) {
            SF_h.push(H_RA); BM_h.push(H_RA * x[i]);
        } else if (x[i] < AD) {
            SF_h.push(H_RA - FtC); BM_h.push(H_RA * x[i] - FtC * (x[i] - AC));
        } else {
            SF_h.push(H_RA - FtC + FtD); BM_h.push(H_RA * x[i] - FtC * (x[i] - AC) + FtD * (x[i] - AD));
        }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'twogear', bearing1: 0, bearing2: AB,
        FtC, FrC, FtD, FrD, AC, AD, AB, V_RA, V_RB, H_RA, H_RB, Wc, Wd,
        details: 'Gear C: Ft=' + FtC.toFixed(2) + 'N, Fr=' + FrC.toFixed(2) + 'N\nGear D: Ft=' + FtD.toFixed(2) + 'N, Fr=' + FrD.toFixed(2) + 'N\nV: RA=' + V_RA.toFixed(2) + 'N, RB=' + V_RB.toFixed(2) + 'N\nH: RA=' + H_RA.toFixed(2) + 'N, RB=' + H_RB.toFixed(2) + 'N'
    };
}

// ========================================
// PULLEY + GEAR SHAFT
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

    const R_pulley = pulleyDia / 2;
    const T2 = T / (R_pulley * (tensionRatio - 1));
    const T1 = tensionRatio * T2;

    const FtB = (2 * T) / gearDia;
    const FrB = FtB * Math.tan(alpha);

    const bearingC = pulleyPos;
    const bearingD = pulleyPos + CD;
    const L = Math.max(bearingD, gearPos) + 50;

    // Vertical: T1+T2+Wa at pulley, FtB+Wb at gear
    const vertPulley = T1 + T2 + Wa;
    const vertGear = FtB + Wb;

    const RDv = (vertPulley * pulleyPos + vertGear * gearPos) / CD;
    const RCv = vertPulley + vertGear - RDv;

    // Horizontal: FrB at gear
    const RDh = (FrB * (gearPos - bearingC)) / CD;
    const RCh = FrB - RDh;

    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        // Vertical
        if (x[i] < bearingC) {
            SF_v.push(-vertPulley); BM_v.push(-vertPulley * x[i]);
        } else if (x[i] < bearingD) {
            SF_v.push(-vertPulley + RCv); BM_v.push(-vertPulley * x[i] + RCv * (x[i] - bearingC));
        } else if (x[i] <= gearPos) {
            SF_v.push(-vertPulley + RCv + RDv); BM_v.push(-vertPulley * x[i] + RCv * (x[i] - bearingC) + RDv * (x[i] - bearingD));
        } else { SF_v.push(0); BM_v.push(0); }

        // Horizontal
        if (x[i] < bearingC) { SF_h.push(0); BM_h.push(0); }
        else if (x[i] < bearingD) { SF_h.push(RCh); BM_h.push(RCh * (x[i] - bearingC)); }
        else if (x[i] <= gearPos) { SF_h.push(RCh + RDh); BM_h.push(RCh * (x[i] - bearingC) + RDh * (x[i] - bearingD)); }
        else { SF_h.push(0); BM_h.push(0); }
    }

    const BM_res = BM_h.map((bh, i) => Math.sqrt(bh * bh + BM_v[i] * BM_v[i]));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'pulleygear', bearing1: bearingC, bearing2: bearingD,
        T1, T2, FtB, FrB, RCv, RDv, RCh, RDh, Wa, Wb, vertPulley, vertGear,
        details: 'T1=' + T1.toFixed(2) + 'N, T2=' + T2.toFixed(2) + 'N\nGear: Ft=' + FtB.toFixed(2) + 'N, Fr=' + FrB.toFixed(2) + 'N\nV: RC=' + RCv.toFixed(2) + 'N, RD=' + RDv.toFixed(2) + 'N\nH: RC=' + RCh.toFixed(2) + 'N, RD=' + RDh.toFixed(2) + 'N'
    };
}

// ========================================
// MULTI PULLEY SHAFT
// ========================================
function calcMultiPulley(T) {
    const AB = parseFloat(document.getElementById('mp_bearingDist').value);
    const posA = parseFloat(document.getElementById('mp_pulleyA_pos').value);
    const diaA = parseFloat(document.getElementById('mp_pulleyA_dia').value);
    const WpA = parseFloat(document.getElementById('mp_pulleyA_weight').value);
    const T1A = parseFloat(document.getElementById('mp_T1A').value);
    const T2A = parseFloat(document.getElementById('mp_T2A').value);

    const posB = parseFloat(document.getElementById('mp_pulleyB_pos').value);
    const diaB = parseFloat(document.getElementById('mp_pulleyB_dia').value);
    const WpB = parseFloat(document.getElementById('mp_pulleyB_weight').value);
    const T1B = parseFloat(document.getElementById('mp_T1B').value);
    const T2B = parseFloat(document.getElementById('mp_T2B').value);

    const forceA = T1A + T2A + WpA;
    const forceB = T1B + T2B + WpB;

    const R2v = (forceA * posA + forceB * posB) / AB;
    const R1v = forceA + forceB - R2v;

    const L = AB;
    const n = 500;
    const x = Array.from({ length: n }, (_, i) => (i / (n - 1)) * L);
    const SF_h = [], BM_h = [], SF_v = [], BM_v = [];

    for (let i = 0; i < n; i++) {
        SF_h.push(0); BM_h.push(0);

        if (x[i] < posA) {
            SF_v.push(R1v); BM_v.push(R1v * x[i]);
        } else if (x[i] < posB) {
            SF_v.push(R1v - forceA); BM_v.push(R1v * x[i] - forceA * (x[i] - posA));
        } else {
            SF_v.push(R1v - forceA - forceB); BM_v.push(R1v * x[i] - forceA * (x[i] - posA) - forceB * (x[i] - posB));
        }
    }

    const BM_res = BM_v.map(b => Math.abs(b));
    const maxBM = Math.max(...BM_res);

    return {
        x, L, SF_h, BM_h, SF_v, BM_v, BM_res, maxBM,
        type: 'multipulley', bearing1: 0, bearing2: AB,
        details: 'Pulley A: T1=' + T1A + 'N, T2=' + T2A + 'N, Force=' + forceA.toFixed(2) + 'N\nPulley B: T1=' + T1B + 'N, T2=' + T2B + 'N, Force=' + forceB.toFixed(2) + 'N\nR1=' + R1v.toFixed(2) + 'N, R2=' + R2v.toFixed(2) + 'N'
    };
}

// ========================================
// DISPLAY RESULTS
// ========================================
function displayResults(T, maxBM, diameter, outerDia, innerDia, shaftType, Cm, Ct, sigmaMax, tauMax, dNormal, dShear, result, hasKeyway) {
    document.getElementById('results').style.display = 'block';
    document.getElementById('diagrams').style.display = 'block';

    document.getElementById('torqueResult').textContent = T.toFixed(2) + ' N-mm (' + (T / 1000).toFixed(2) + ' N-m)';
    document.getElementById('momentResult').textContent = maxBM.toFixed(2) + ' N-mm';

    if (shaftType === 'solid') {
        document.getElementById('diameterResult').textContent = 'd = ' + outerDia + ' mm';
    } else {
        document.getElementById('diameterResult').textContent = 'do = ' + outerDia + ' mm, di = ' + innerDia.toFixed(2) + ' mm';
    }

    let html = '<h3>📋 Step-by-Step Solution</h3>';
    html += '<p><strong>Torque:</strong> T = (9.55 × 10⁶ × ' + document.getElementById('power').value + ') / ' + document.getElementById('speed').value + ' = ' + T.toFixed(2) + ' N-mm</p>';
    html += '<p><strong>σ_max = ' + sigmaMax.toFixed(2) + ' MPa</strong></p>';
    html += '<p><strong>τ_max = ' + tauMax.toFixed(2) + ' MPa</strong></p>';
    if (hasKeyway) html += '<p><strong>⚠ Keyway reduction applied (× 0.75)</strong></p>';
    html += '<p><strong>Cm = ' + Cm.toFixed(2) + ', Ct = ' + Ct.toFixed(2) + '</strong></p>';

    html += '<h3>📐 Force Analysis</h3>';
    html += '<pre style="background:#f0f0f0; padding:10px; border-radius:5px; white-space:pre-wrap;">' + result.details + '</pre>';

    html += '<h3>📊 Bending Moment</h3>';
    html += '<p><strong>Maximum Bending Moment M = ' + maxBM.toFixed(2) + ' N-mm</strong></p>';

    html += '<h3>📏 Diameter Calculation</h3>';
    html += '<p><strong>Eq.(i) Max Normal Stress Theory:</strong></p>';
    html += '<p>d = [16/(π×' + sigmaMax.toFixed(2) + ') × {' + Cm + '×' + maxBM.toFixed(0) + ' + √((' + Cm + '×' + maxBM.toFixed(0) + ')² + (' + Ct + '×' + T.toFixed(0) + ')²)}]^(1/3)</p>';
    html += '<p><strong>d = ' + dNormal.toFixed(2) + ' mm ...Eq.(i)</strong></p>';

    html += '<p><strong>Eq.(ii) Max Shear Stress Theory:</strong></p>';
    html += '<p>d = [16/(π×' + tauMax.toFixed(2) + ') × √((' + Cm + '×' + maxBM.toFixed(0) + ')² + (' + Ct + '×' + T.toFixed(0) + ')²)]^(1/3)</p>';
    html += '<p><strong>d = ' + dShear.toFixed(2) + ' mm ...Eq.(ii)</strong></p>';

    html += '<p>Select maximum: d = ' + Math.max(dNormal, dShear).toFixed(2) + ' mm</p>';
    html += '<p style="color:red; font-size:1.2em; font-weight:bold;">∴ Standard size of shaft, d = ' + outerDia + ' mm</p>';

    document.getElementById('detailedResults').innerHTML = html;
    document.getElementById('results').scrollIntoView({ behavior: 'smooth' });
}

// ========================================
// GENERATE DIAGRAMS
// ========================================
function generateDiagrams(result, T) {
    const x = result.x;
    const n = x.length;
    const torqueData = Array(n).fill(T);

    document.getElementById('diagrams').innerHTML =
        '<h2>📈 Horizontal Loading Diagram (HLD)</h2>' +
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
    const idxH = result.BM_h.findIndex(b => Math.abs(b) === maxH);
    Plotly.newPlot('bmd_h', [{ x: x, y: result.BM_h, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(219,64,82)', width: 3 } }],
        { title: 'Bending Moment Diagram - Horizontal Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'BM (N-mm)' },
            annotations: [{ x: x[idxH], y: result.BM_h[idxH], text: 'M_max = ' + result.BM_h[idxH].toFixed(0) + ' N-mm', showarrow: true, arrowhead: 2, ax: 40, ay: -40 }] });

    // VLD - SFD
    Plotly.newPlot('sfd_v', [{ x: x, y: result.SF_v, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(55,128,191)', width: 3 } }],
        { title: 'Shear Force Diagram - Vertical Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'SF (N)' } });

    // VLD - BMD
    const maxV = Math.max(...result.BM_v.map(Math.abs));
    const idxV = result.BM_v.findIndex(b => Math.abs(b) === maxV);
    Plotly.newPlot('bmd_v', [{ x: x, y: result.BM_v, type: 'scatter', fill: 'tozeroy', line: { color: 'rgb(219,64,82)', width: 3 } }],
        { title: 'Bending Moment Diagram - Vertical Plane', xaxis: { title: 'Distance (mm)' }, yaxis: { title: 'BM (N-mm)
