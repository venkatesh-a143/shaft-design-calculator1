function toggleProblemType(){
    var v=document.querySelector('input[name="problemType"]:checked').value;
    document.getElementById('normalSection').style.display=v==='normal'?'block':'none';
    document.getElementById('udlSection').style.display=v==='udl'?'block':'none';
    document.getElementById('singlePulleySection').style.display=v==='singlepulley'?'block':'none';
    document.getElementById('twoGearSection').style.display=v==='twogear'?'block':'none';
    document.getElementById('pulleyGearSection').style.display=v==='pulleygear'?'block':'none';
    document.getElementById('multiPulleySection').style.display=v==='multipulley'?'block':'none';
}
function toggleStressMethod(){
    var m=document.querySelector('input[name="stressMethod"]:checked').value;
    document.getElementById('directStressInputs').style.display=m==='direct'?'block':'none';
    document.getElementById('materialInputs').style.display=m==='material'?'block':'none';
    document.getElementById('workingStressInputs').style.display=m==='working'?'block':'none';
    document.getElementById('yieldStressInputs').style.display=m==='yield'?'block':'none';
}
function toggleShaftType(){document.getElementById('hollowInputs').style.display=document.querySelector('input[name="shaftType"]:checked').value==='hollow'?'block':'none';}
function toggleTensionInput(){
    var v=document.querySelector('input[name="tensionInput"]:checked').value;
    document.getElementById('directTensionInputs').style.display=v==='direct'?'block':'none';
    document.getElementById('ratioTensionInputs').style.display=v==='ratio'?'block':'none';
}
function toggleModulusG(){document.getElementById('modulusGInput').style.display=document.getElementById('hasModulusG').checked?'block':'none';}
function toggleShaftLength(){document.getElementById('shaftLengthInput').style.display=document.getElementById('hasShaftLength').checked?'block':'none';}
function toggleAngleTwist(){document.getElementById('angleTwistInput').style.display=document.getElementById('hasAngleTwist').checked?'block':'none';}
document.addEventListener('DOMContentLoaded',function(){
    document.getElementById('calcDeflection').addEventListener('change',function(){
        document.getElementById('deflectionInputs').style.display=this.checked?'block':'none';
    });
});
function evalFraction(s){s=s.trim();if(s.indexOf('/')!==-1){var p=s.split('/');return parseFloat(p[0])/parseFloat(p[1]);}return parseFloat(s);}

function calculateShaft(){
    var problemType=document.querySelector('input[name="problemType"]:checked').value;
    var P=parseFloat(document.getElementById('power').value);
    var N=parseFloat(document.getElementById('speed').value);
    var T=(9.55e6*P)/N;

    var stressMethod=document.querySelector('input[name="stressMethod"]:checked').value;
    var sigmaMax,tauMax;
    if(stressMethod==='direct'){tauMax=parseFloat(document.getElementById('tauAllow').value);sigmaMax=2*tauMax;}
    else if(stressMethod==='material'){
        var sigUt=parseFloat(document.getElementById('sigmaUt').value);var sigYt=parseFloat(document.getElementById('sigmaYt').value);var nP=parseFloat(document.getElementById('nPrime').value);
        var sigU=sigUt/nP;var sigY=sigYt/nP;sigmaMax=Math.min(0.36*sigU,0.6*sigY);tauMax=Math.min(0.18*sigU,0.3*sigY);
    }else if(stressMethod==='working'){sigmaMax=parseFloat(document.getElementById('sigmaWorking').value);tauMax=parseFloat(document.getElementById('tauWorking').value);}
    else if(stressMethod==='yield'){var sYt=parseFloat(document.getElementById('sigmaYtYield').value);var fos=parseFloat(document.getElementById('fosYield').value);sigmaMax=sYt/fos;tauMax=sigmaMax/2;}

    if(document.getElementById('hasKeyway').checked){sigmaMax=0.75*sigmaMax;tauMax=0.75*tauMax;}
    var Cm=parseFloat(document.getElementById('cm').value);var Ct=parseFloat(document.getElementById('ct').value);

    var result;
    if(problemType==='normal') result=calcNormalShaft(T);
    else if(problemType==='udl') result=calcUDL(T);
    else if(problemType==='singlepulley') result=calcSinglePulley(T);
    else if(problemType==='twogear') result=calcTwoGear(T);
    else if(problemType==='pulleygear') result=calcPulleyGear(T);
    else if(problemType==='multipulley') result=calcMultiPulley(T);

    var shaftType=document.querySelector('input[name="shaftType"]:checked').value;
    var maxBM=result.maxBM;
    var term1=Cm*maxBM;var term2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));
    var dNormal=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2),1/3);
    var dShear=Math.pow((16/(Math.PI*tauMax))*term2,1/3);

    var outerDia,innerDia,diameter;
    var stdSizes=[10,12,16,20,25,30,35,40,45,50,55,56,60,63,65,70,75,80,85,90,95,100,110,120];
    if(shaftType==='solid'){
        diameter=Math.max(dNormal,dShear);outerDia=stdSizes.find(function(s){return s>=diameter;})||Math.ceil(diameter/10)*10;innerDia=0;
    }else{
        var k=evalFraction(document.getElementById('diameterRatio').value);
        var doNorm=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2)*(1/(1-Math.pow(k,4))),1/3);
        var doShear2=Math.pow((16/(Math.PI*tauMax))*term2*(1/(1-Math.pow(k,4))),1/3);
        outerDia=Math.max(doNorm,doShear2);outerDia=stdSizes.find(function(s){return s>=outerDia;})||Math.ceil(outerDia/10)*10;
        innerDia=k*outerDia;diameter=outerDia;
    }

    var deflectionText='';
    if(document.getElementById('calcDeflection').checked){
        var Gd=parseFloat(document.getElementById('deflG').value)*1000;var Ld=parseFloat(document.getElementById('deflL').value);
        var D4=Math.pow(outerDia,4);var theta=(584*T*Ld)/(Gd*D4);
        deflectionText='<p>\u03B8 = (584 \u00D7 T \u00D7 L) / (G \u00D7 d\u2074)</p><p>= (584 \u00D7 '+T.toFixed(2)+' \u00D7 '+Ld+') / ('+Gd+' \u00D7 '+outerDia+'\u2074)</p><p><strong>\u03B8 = '+theta.toFixed(4)+' degrees</strong></p>';
    }

    displayResults(T,maxBM,diameter,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflectionText);
    generateDiagrams(result,T);
}

function calcNormalShaft(T){
    var AB=parseFloat(document.getElementById('normal_bearingDist').value);
    var gp=parseFloat(document.getElementById('normal_gearPos').value);
    var gd=parseFloat(document.getElementById('normal_gearDia').value);
    var alpha=parseFloat(document.getElementById('normal_pressureAngle').value)*Math.PI/180;
    var Wg=parseFloat(document.getElementById('normal_gearWeight').value);
    var Ft=(2*T)/gd;var Fr=Ft*Math.tan(alpha);
    var RAh=Ft*(AB-gp)/AB;var RBh=Ft-RAh;
    var vf=Fr+Wg;var RAv=vf*(AB-gp)/AB;var RBv=vf-RAv;
    var L=AB;var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);
        if(xi<gp){SF_h.push(RAh);BM_h.push(RAh*xi);SF_v.push(RAv);BM_v.push(RAv*xi);}
        else{SF_h.push(RAh-Ft);BM_h.push(RAh*xi-Ft*(xi-gp));SF_v.push(RAv-vf);BM_v.push(RAv*xi-vf*(xi-gp));}}
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'normal',
        details:'F<sub>t</sub> = 2T/D = '+Ft.toFixed(2)+' N\nF<sub>r</sub> = F<sub>t</sub> \u00D7 tan(\u03B1) = '+Fr.toFixed(2)+' N\n\nHorizontal Reactions:\nR<sub>AH</sub> = '+RAh.toFixed(2)+' N\nR<sub>BH</sub> = '+RBh.toFixed(2)+' N\n\nVertical Reactions:\nR<sub>AV</sub> = '+RAv.toFixed(2)+' N\nR<sub>BV</sub> = '+RBv.toFixed(2)+' N'};
}

function calcUDL(T){
    var L=parseFloat(document.getElementById('totalLength').value);var lU=parseFloat(document.getElementById('udlLength').value);
    var w=parseFloat(document.getElementById('udlIntensity').value);var b1=(L-lU)/2;var b2=b1+lU;var RA=(w*lU)/2;
    var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);
        if(xi<b1||xi>b2){SF_h.push(0);BM_h.push(0);SF_v.push(0);BM_v.push(0);}
        else{var d=xi-b1;var sf=RA-w*d;var bm=RA*d-w*d*d/2;SF_h.push(sf);BM_h.push(bm);SF_v.push(sf);BM_v.push(bm);}}
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'udl',
        details:'w = '+w+' N/mm, l = '+lU+' mm\nTotal Load W = w \u00D7 l = '+(w*lU).toFixed(2)+' N\nR<sub>A</sub> = R<sub>B</sub> = W/2 = '+(w*lU/2).toFixed(2)+' N\nM<sub>max</sub> at centre = wl\u00B2/8 = '+(w*lU*lU/8).toFixed(2)+' N\u00B7mm'};
}

function calcSinglePulley(T){
    var AB=parseFloat(document.getElementById('sp_bearingDist').value);
    var pp=parseFloat(document.getElementById('sp_pulleyPos').value);
    var pd=parseFloat(document.getElementById('sp_pulleyDia').value);
    var Wp=parseFloat(document.getElementById('sp_pulleyWeight').value);var T1,T2;
    if(document.querySelector('input[name="tensionInput"]:checked').value==='direct'){T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);}
    else{var ratio=parseFloat(document.getElementById('sp_tensionRatio').value);T2=T/((pd/2)*(ratio-1));T1=ratio*T2;}
    var bf=T1+T2+Wp;var RAv=bf*(AB-pp)/AB;var RBv=bf-RAv;
    var L=AB;var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);SF_h.push(0);BM_h.push(0);
        if(xi<pp){SF_v.push(RAv);BM_v.push(RAv*xi);}
        else{SF_v.push(RAv-bf);BM_v.push(RAv*xi-bf*(xi-pp));}}
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.abs(BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'singlepulley',
        details:'T\u2081 = '+T1.toFixed(2)+' N\nT\u2082 = '+T2.toFixed(2)+' N\nTotal Belt Force (T\u2081+T\u2082) = '+(T1+T2).toFixed(2)+' N\nPulley Weight W = '+Wp+' N\nTotal Vertical Force = '+bf.toFixed(2)+' N\n\nR<sub>A</sub> = '+RAv.toFixed(2)+' N\nR<sub>B</sub> = '+RBv.toFixed(2)+' N'};
}

// ========================================
// TWO GEARS - FIXED R_BV FORMULA
// ========================================
// From textbook: Taking moments about A, equating to zero:
// R_BV x AB = F_tD(AD) - F_rC(AC)
// R_BV = [F_tD x AD - F_rC x AC] / AB
// R_AV + R_BV = F_tD - F_rC
// R_AV = F_tD - F_rC - R_BV
function calcTwoGear(T){
    var AB=parseFloat(document.getElementById('bearingDistance').value);
    var AC=parseFloat(document.getElementById('gearC_pos').value);
    var AD=parseFloat(document.getElementById('gearD_pos').value);
    var Dc=parseFloat(document.getElementById('gearC_dia').value);
    var Dd=parseFloat(document.getElementById('gearD_dia').value);
    var Wc=parseFloat(document.getElementById('gearC_weight').value);
    var Wd=parseFloat(document.getElementById('gearD_weight').value);
    var alpha=parseFloat(document.getElementById('twoGear_pressureAngle').value)*Math.PI/180;

    var FtC=(2*T)/Dc;var FrC=FtC*Math.tan(alpha);
    var FtD=(2*T)/Dd;var FrD=FtD*Math.tan(alpha);

    // ===== VERTICAL PLANE (Radial forces) =====
    // Taking moments about bearing A and equating to zero:
    // R_BV x AB = F_tD(AD) - F_rC(AC)
    // (Using textbook formula exactly as shown in images)
    var V_RB = (FtD*AD - FrC*AC) / AB;
    // Also: R_AV + R_BV = F_tD - F_rC
    var V_RA = FtD - FrC - V_RB;

    // If weights are given, add them
    if(Wc > 0 || Wd > 0){
        V_RB = ((FtD+Wd)*AD - (FrC-Wc)*AC) / AB;
        V_RA = (FtD+Wd) - (FrC-Wc) - V_RB;
    }

    // ===== HORIZONTAL PLANE (Tangential forces) =====
    // Taking moments about bearing A:
    // R_BH x AB + F_tD(AD) = F_tC(AC)
    // R_BH = [F_tC x AC - F_tD x AD] / AB
    var H_RB = (FtC*AC - FtD*AD) / AB;
    // Check if negative (direction reversal as per textbook)
    var H_RB_negative = H_RB < 0;
    // R_AH + F_tD = F_tC + R_BH (from revised HLD)
    var H_RA = FtC + H_RB - FtD;
    // If R_AH is negative, direction reverses
    var H_RA_negative = H_RA < 0;

    // Moments at C
    var M_CV = V_RA * AC;
    var M_CH = H_RA * AC;

    // Moments at D
    var M_DV = V_RA * AD - (Wc > 0 ? (FrC-Wc) : FrC) * (AD - AC);
    // For horizontal: need to use correct formula based on sign
    var M_DH;
    if(H_RB_negative){
        M_DH = -Math.abs(H_RB) * (AB - AD);
    }else{
        M_DH = H_RA * AD - FtC * (AD - AC);
    }

    // Resultant moments
    var M_C = Math.sqrt(M_CV*M_CV + M_CH*M_CH);
    var M_D = Math.sqrt(M_DV*M_DV + M_DH*M_DH);

    var L=AB;var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);
        // Vertical plane
        if(xi<AC){SF_v.push(V_RA);BM_v.push(V_RA*xi);}
        else if(xi<AD){
            var vForceC = Wc > 0 ? (FrC - Wc) : FrC;
            SF_v.push(V_RA - vForceC);BM_v.push(V_RA*xi - vForceC*(xi-AC));
        }
        else{
            var vForceC2 = Wc > 0 ? (FrC - Wc) : FrC;
            var vForceD = Wd > 0 ? (FtD + Wd) : FtD;
            SF_v.push(V_RA - vForceC2 + vForceD);BM_v.push(V_RA*xi - vForceC2*(xi-AC) + vForceD*(xi-AD));
        }
        // Horizontal plane
        if(xi<AC){SF_h.push(H_RA);BM_h.push(H_RA*xi);}
        else if(xi<AD){SF_h.push(H_RA-FtC);BM_h.push(H_RA*xi-FtC*(xi-AC));}
        else{SF_h.push(H_RA-FtC+FtD);BM_h.push(H_RA*xi-FtC*(xi-AC)+FtD*(xi-AD));}
    }

    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);

    var details = '';
    details += 'Gear C: F<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = '+FtC.toFixed(2)+' N\n';
    details += 'F<sub>rC</sub> = F<sub>tC</sub> \u00D7 tan \u03B1 = '+FrC.toFixed(2)+' N\n\n';
    details += 'Gear D: F<sub>tD</sub> = 2M<sub>t</sub>/D<sub>D</sub> = '+FtD.toFixed(2)+' N\n';
    details += 'F<sub>rD</sub> = F<sub>tD</sub> \u00D7 tan \u03B1 = '+FrD.toFixed(2)+' N\n\n';
    details += '\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\n';
    details += 'Taking moments about A:\n';
    details += 'R<sub>BV</sub> \u00D7 '+AB+' = F<sub>tD</sub>(AD) \u2212 F<sub>rC</sub>(AC)\n';
    details += 'R<sub>BV</sub> \u00D7 '+AB+' = ('+FtD.toFixed(2)+' \u00D7 '+AD+') \u2212 ('+FrC.toFixed(2)+' \u00D7 '+AC+')\n';
    details += 'R<sub>BV</sub> = '+V_RB.toFixed(2)+' N\n\n';
    details += 'Also R<sub>AV</sub> + R<sub>BV</sub> = F<sub>tD</sub> \u2212 F<sub>rC</sub>\n';
    details += 'R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';
    details += 'Moments:\n';
    details += 'M<sub>AV</sub> = M<sub>BV</sub> = 0\n';
    details += 'M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    details += 'M<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+(V_RB*(AB-AD)).toFixed(2)+' N\u00B7mm\n\n';
    details += '\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\n';
    details += 'Taking moments about A:\n';
    details += 'R<sub>BH</sub> \u00D7 '+AB+' + F<sub>tD</sub>(AD) = F<sub>tC</sub>(AC)\n';
    details += 'R<sub>BH</sub> = '+H_RB.toFixed(2)+' N';
    if(H_RB_negative) details += ' (acts in opposite direction)';
    details += '\n\nR<sub>AH</sub> = '+H_RA.toFixed(2)+' N';
    if(H_RA_negative) details += ' (acts in opposite direction)';
    details += '\n\nMoments:\n';
    details += 'M<sub>AH</sub> = M<sub>BH</sub> = 0\n';
    details += 'M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    details += 'M<sub>DH</sub> = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';
    details += '\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\n';
    details += 'M<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + '+M_CH.toFixed(2)+'\u00B2)\n';
    details += 'M<sub>C</sub> = '+M_C.toFixed(2)+' N\u00B7mm\n\n';
    details += 'M<sub>D</sub> = \u221A(M\u00B2<sub>DV</sub> + M\u00B2<sub>DH</sub>) = \u221A('+M_DV.toFixed(2)+'\u00B2 + '+M_DH.toFixed(2)+'\u00B2)\n';
    details += 'M<sub>D</sub> = '+M_D.toFixed(2)+' N\u00B7mm\n\n';
    details += 'Maximum BM occurs at '+(M_C>M_D?'C':'D')+'\n';
    details += '<strong>M = M<sub>'+(M_C>M_D?'C':'D')+'</sub> = '+Math.max(M_C,M_D).toFixed(2)+' N\u00B7mm</strong>';

    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'twogear',details:details};
}

function calcPulleyGear(T){
    var CD=parseFloat(document.getElementById('pg_bearingDistance').value);
    var pp=parseFloat(document.getElementById('pulley_pos').value);
    var pd=parseFloat(document.getElementById('pulley_dia').value);
    var Wa=parseFloat(document.getElementById('pulley_weight').value);
    var tr=parseFloat(document.getElementById('tensionRatio').value);
    var gp=parseFloat(document.getElementById('gear_pos').value);
    var gd=parseFloat(document.getElementById('gear_dia').value);
    var Wb=parseFloat(document.getElementById('gear_weight').value);
    var alpha=parseFloat(document.getElementById('pg_pressureAngle').value)*Math.PI/180;

    var Rp=pd/2;var T2=T/(Rp*(tr-1));var T1=tr*T2;
    var FtB=(2*T)/gd;var FrB=FtB*Math.tan(alpha);

    var bC=pp;var bD=pp+CD;var L=Math.max(bD,gp)+50;
    var vP=T1+T2+Wa;var vG=FtB+Wb;
    var RDv=(vP*pp+vG*gp)/CD;var RCv=vP+vG-RDv;
    var RDh=(FrB*(gp-bC))/CD;var RCh=FrB-RDh;

    var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);
        if(xi<bC){SF_v.push(-vP);BM_v.push(-vP*xi);SF_h.push(0);BM_h.push(0);}
        else if(xi<bD){SF_v.push(-vP+RCv);BM_v.push(-vP*xi+RCv*(xi-bC));SF_h.push(RCh);BM_h.push(RCh*(xi-bC));}
        else if(xi<=gp){SF_v.push(-vP+RCv+RDv);BM_v.push(-vP*xi+RCv*(xi-bC)+RDv*(xi-bD));SF_h.push(RCh+RDh);BM_h.push(RCh*(xi-bC)+RDh*(xi-bD));}
        else{SF_v.push(0);BM_v.push(0);SF_h.push(0);BM_h.push(0);}}
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);

    var details = '';
    details += 'T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';
    details += 'Gear: F<sub>t</sub> = '+FtB.toFixed(2)+' N, F<sub>r</sub> = '+FrB.toFixed(2)+' N\n\n';
    details += 'Vertical Reactions:\nR<sub>CV</sub> = '+RCv.toFixed(2)+' N\nR<sub>DV</sub> = '+RDv.toFixed(2)+' N\n\n';
    details += 'Horizontal Reactions:\nR<sub>CH</sub> = '+RCh.toFixed(2)+' N\nR<sub>DH</sub> = '+RDh.toFixed(2)+' N';
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'pulleygear',details:details};
}

function calcMultiPulley(T){
    var AB=parseFloat(document.getElementById('mp_bearingDist').value);
    var pA=parseFloat(document.getElementById('mp_pulleyA_pos').value);var dA=parseFloat(document.getElementById('mp_pulleyA_dia').value);
    var WpA=parseFloat(document.getElementById('mp_pulleyA_weight').value);
    var T1A=parseFloat(document.getElementById('mp_T1A').value);var T2A=parseFloat(document.getElementById('mp_T2A').value);
    var pB=parseFloat(document.getElementById('mp_pulleyB_pos').value);var dB=parseFloat(document.getElementById('mp_pulleyB_dia').value);
    var WpB=parseFloat(document.getElementById('mp_pulleyB_weight').value);
    var T1B=parseFloat(document.getElementById('mp_T1B').value);var T2B=parseFloat(document.getElementById('mp_T2B').value);

    var fA=T1A+T2A+WpA;var fB=T1B+T2B+WpB;
    var R2v=(fA*pA+fB*pB)/AB;var R1v=fA+fB-R2v;
    var L=AB;var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);SF_h.push(0);BM_h.push(0);
        if(xi<pA){SF_v.push(R1v);BM_v.push(R1v*xi);}
        else if(xi<pB){SF_v.push(R1v-fA);BM_v.push(R1v*xi-fA*(xi-pA));}
        else{SF_v.push(R1v-fA-fB);BM_v.push(R1v*xi-fA*(xi-pA)-fB*(xi-pB));}}
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.abs(BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);

    var details = '';
    details += 'Pulley A: T\u2081 = '+T1A+' N, T\u2082 = '+T2A+' N, W = '+WpA+' N\n';
    details += 'Total Force at A = '+fA.toFixed(2)+' N\n\n';
    details += 'Pulley B: T\u2081 = '+T1B+' N, T\u2082 = '+T2B+' N, W = '+WpB+' N\n';
    details += 'Total Force at B = '+fB.toFixed(2)+' N\n\n';
    details += 'R\u2081 = '+R1v.toFixed(2)+' N\nR\u2082 = '+R2v.toFixed(2)+' N';
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'multipulley',details:details};
}

function displayResults(T,maxBM,diameter,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflectionText){
    document.getElementById('results').style.display='block';
    document.getElementById('diagrams').style.display='block';
    document.getElementById('torqueResult').textContent=T.toFixed(2)+' N\u00B7mm ('+(T/1000).toFixed(2)+' N\u00B7m)';
    document.getElementById('momentResult').textContent=maxBM.toFixed(2)+' N\u00B7mm';
    if(shaftType==='solid') document.getElementById('diameterResult').textContent='d = '+outerDia+' mm';
    else document.getElementById('diameterResult').textContent='d\u2092 = '+outerDia+' mm, d\u1D62 = '+innerDia.toFixed(2)+' mm';

    var h='<h3>Step-by-Step Solution</h3>';
    h+='<p><strong>Torque:</strong> T = (9.55 \u00D7 10\u2076 \u00D7 P) / n = (9.55 \u00D7 10\u2076 \u00D7 '+document.getElementById('power').value+') / '+document.getElementById('speed').value+' = '+T.toFixed(2)+' N\u00B7mm</p>';
    h+='<p><strong>\u03C3<sub>max</sub> = '+sigmaMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>\u03C4<sub>max</sub> = '+tauMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>C<sub>m</sub> = '+Cm.toFixed(2)+', C<sub>t</sub> = '+Ct.toFixed(2)+'</strong></p>';
    h+='<h3>Force Analysis</h3>';
    h+='<pre style="background:#f0f0f0;padding:15px;border-radius:5px;white-space:pre-wrap;font-size:14px;line-height:1.8;">'+result.details+'</pre>';
    h+='<h3>Maximum Bending Moment</h3>';
    h+='<p><strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong></p>';
    h+='<h3>Diameter Calculation</h3>';
    h+='<p><strong>Eq.(i) Max Normal Stress Theory (3.6(a)/Pg 51, DHB):</strong></p>';
    h+='<p>d = [16/(\u03C0\u03C3<sub>max</sub>) \u00D7 {C<sub>m</sub>M + \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)}]<sup>1/3</sup></p>';
    h+='<p>d = [16/(\u03C0 \u00D7 '+sigmaMax.toFixed(2)+') \u00D7 {'+Cm+' \u00D7 '+maxBM.toFixed(0)+' + \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(0)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(0)+')\u00B2)}]<sup>1/3</sup></p>';
    h+='<p><strong>d = '+dNormal.toFixed(2)+' mm ...Eq.(i)</strong></p>';
    h+='<p><strong>Eq.(ii) Max Shear Stress Theory (3.6(b)/Pg 51, DHB):</strong></p>';
    h+='<p>d = [16/(\u03C0\u03C4<sub>max</sub>) \u00D7 \u221A((C<sub>m</sub>M)\u00B2 + (C<sub
