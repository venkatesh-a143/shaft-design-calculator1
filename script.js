// ========== TOGGLE FUNCTIONS ==========
function toggleProblemType(){
    var v=document.querySelector('input[name="problemType"]:checked').value;
    var sections=['normalSection','udlSection','singlePulleySection','twoGearSection','pulleyGearSection','multiPulleySection'];
    var map={normal:'normalSection',udl:'udlSection',singlepulley:'singlePulleySection',twogear:'twoGearSection',pulleygear:'pulleyGearSection',multipulley:'multiPulleySection'};
    for(var i=0;i<sections.length;i++){document.getElementById(sections[i]).style.display='none';}
    document.getElementById(map[v]).style.display='block';
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

// ========== MAIN CALCULATION ==========
function calculateShaft(){
    var problemType=document.querySelector('input[name="problemType"]:checked').value;
    var P=parseFloat(document.getElementById('power').value);
    var N=parseFloat(document.getElementById('speed').value);
    var T=(9.55e6*P)/N;

    var stressMethod=document.querySelector('input[name="stressMethod"]:checked').value;
    var sigmaMax,tauMax;
    if(stressMethod==='direct'){tauMax=parseFloat(document.getElementById('tauAllow').value);sigmaMax=2*tauMax;}
    else if(stressMethod==='material'){
        var sigUt=parseFloat(document.getElementById('sigmaUt').value);var sigYt=parseFloat(document.getElementById('sigmaYt').value);
        var nP=parseFloat(document.getElementById('nPrime').value);var sigU=sigUt/nP;var sigY=sigYt/nP;
        sigmaMax=Math.min(0.36*sigU,0.6*sigY);tauMax=Math.min(0.18*sigU,0.3*sigY);
    }else if(stressMethod==='working'){
        sigmaMax=parseFloat(document.getElementById('sigmaWorking').value);tauMax=parseFloat(document.getElementById('tauWorking').value);
    }else if(stressMethod==='yield'){
        var sYt=parseFloat(document.getElementById('sigmaYtYield').value);var fos=parseFloat(document.getElementById('fosYield').value);
        sigmaMax=sYt/fos;tauMax=sigmaMax/2;
    }
    if(document.getElementById('hasKeyway').checked){sigmaMax=0.75*sigmaMax;tauMax=0.75*tauMax;}
    var Cm=parseFloat(document.getElementById('cm').value);
    var Ct=parseFloat(document.getElementById('ct').value);

    var result;
    if(problemType==='normal')result=calcNormalShaft(T);
    else if(problemType==='udl')result=calcUDL(T);
    else if(problemType==='singlepulley')result=calcSinglePulley(T);
    else if(problemType==='twogear')result=calcTwoGear(T);
    else if(problemType==='pulleygear')result=calcPulleyGear(T);
    else if(problemType==='multipulley')result=calcMultiPulley(T);

    var shaftType=document.querySelector('input[name="shaftType"]:checked').value;
    var maxBM=result.maxBM;
    var term1=Cm*maxBM;
    var term2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));
    var dNormal=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2),1/3);
    var dShear=Math.pow((16/(Math.PI*tauMax))*term2,1/3);

    var outerDia,innerDia;
    var stdSizes=[10,12,16,20,25,30,35,40,45,50,55,56,60,63,65,70,75,80,85,90,95,100,110,120];
    if(shaftType==='solid'){
        var d=Math.max(dNormal,dShear);
        outerDia=stdSizes.find(function(s){return s>=d;})||Math.ceil(d/10)*10;
        innerDia=0;
    }else{
        var k=evalFraction(document.getElementById('diameterRatio').value);
        var doN=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2)*(1/(1-Math.pow(k,4))),1/3);
        var doS=Math.pow((16/(Math.PI*tauMax))*term2*(1/(1-Math.pow(k,4))),1/3);
        outerDia=Math.max(doN,doS);
        outerDia=stdSizes.find(function(s){return s>=outerDia;})||Math.ceil(outerDia/10)*10;
        innerDia=k*outerDia;
    }

    var deflText='';
    if(document.getElementById('calcDeflection').checked){
        var Gd=parseFloat(document.getElementById('deflG').value)*1000;
        var Ld=parseFloat(document.getElementById('deflL').value);
        var D4=Math.pow(outerDia,4);
        var theta=(584*T*Ld)/(Gd*D4);
        deflText='<p>\u03B8 = (584 \u00D7 T \u00D7 L) / (G \u00D7 d\u2074)</p>';
        deflText+='<p>= (584 \u00D7 '+T.toFixed(2)+' \u00D7 '+Ld+') / ('+Gd+' \u00D7 '+outerDia+'\u2074)</p>';
        deflText+='<p><strong>\u03B8 = '+theta.toFixed(4)+' degrees</strong></p>';
    }

    displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText);
    generateDiagrams(result,T);
}

// ========== NORMAL SHAFT ==========
function calcNormalShaft(T){
    var AB=parseFloat(document.getElementById('normal_bearingDist').value);
    var gp=parseFloat(document.getElementById('normal_gearPos').value);
    var gd=parseFloat(document.getElementById('normal_gearDia').value);
    var alpha=parseFloat(document.getElementById('normal_pressureAngle').value)*Math.PI/180;
    var Wg=parseFloat(document.getElementById('normal_gearWeight').value);
    var Ft=(2*T)/gd;var Fr=Ft*Math.tan(alpha);
    var RAh=Ft*(AB-gp)/AB;var RBh=Ft-RAh;
    var vf=Fr+Wg;var RAv=vf*(AB-gp)/AB;var RBv=vf-RAv;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);
        if(xi<gp){SFh.push(RAh);BMh.push(RAh*xi);SFv.push(RAv);BMv.push(RAv*xi);}
        else{SFh.push(RAh-Ft);BMh.push(RAh*xi-Ft*(xi-gp));SFv.push(RAv-vf);BMv.push(RAv*xi-vf*(xi-gp));}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:Math.max.apply(null,BR),
        details:'F<sub>t</sub> = 2T/D = 2 \u00D7 '+T.toFixed(2)+' / '+gd+' = '+Ft.toFixed(2)+' N\nF<sub>r</sub> = F<sub>t</sub> \u00D7 tan \u03B1 = '+Ft.toFixed(2)+' \u00D7 tan '+document.getElementById('normal_pressureAngle').value+'\u00B0 = '+Fr.toFixed(2)+' N\n\nHorizontal Reactions:\nR<sub>AH</sub> = '+RAh.toFixed(2)+' N\nR<sub>BH</sub> = '+RBh.toFixed(2)+' N\n\nVertical Reactions:\nR<sub>AV</sub> = '+RAv.toFixed(2)+' N\nR<sub>BV</sub> = '+RBv.toFixed(2)+' N'};
}

// ========== UDL SHAFT ==========
function calcUDL(T){
    var L=parseFloat(document.getElementById('totalLength').value);
    var lU=parseFloat(document.getElementById('udlLength').value);
    var w=parseFloat(document.getElementById('udlIntensity').value);
    var b1=(L-lU)/2;var b2=b1+lU;var RA=(w*lU)/2;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);
        if(xi<b1||xi>b2){SFh.push(0);BMh.push(0);SFv.push(0);BMv.push(0);}
        else{var d=xi-b1;SFh.push(RA-w*d);BMh.push(RA*d-w*d*d/2);SFv.push(RA-w*d);BMv.push(RA*d-w*d*d/2);}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:Math.max.apply(null,BR),
        details:'w = '+w+' N/mm, l = '+lU+' mm\nW = w \u00D7 l = '+(w*lU).toFixed(2)+' N\nR<sub>A</sub> = R<sub>B</sub> = W/2 = '+RA.toFixed(2)+' N\nM<sub>max</sub> = wl\u00B2/8 = '+(w*lU*lU/8).toFixed(2)+' N\u00B7mm'};
}

// ========== SINGLE PULLEY ==========
function calcSinglePulley(T){
    var AB=parseFloat(document.getElementById('sp_bearingDist').value);
    var pp=parseFloat(document.getElementById('sp_pulleyPos').value);
    var pd=parseFloat(document.getElementById('sp_pulleyDia').value);
    var Wp=parseFloat(document.getElementById('sp_pulleyWeight').value);
    var T1,T2;
    if(document.querySelector('input[name="tensionInput"]:checked').value==='direct'){
        T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);
    }else{var r=parseFloat(document.getElementById('sp_tensionRatio').value);T2=T/((pd/2)*(r-1));T1=r*T2;}
    var bf=T1+T2+Wp;var RAv=bf*(AB-pp)/AB;var RBv=bf-RAv;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);SFh.push(0);BMh.push(0);
        if(xi<pp){SFv.push(RAv);BMv.push(RAv*xi);}
        else{SFv.push(RAv-bf);BMv.push(RAv*xi-bf*(xi-pp));}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.abs(BMv[i]));
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:Math.max.apply(null,BR),
        details:'T\u2081 = '+T1.toFixed(2)+' N\nT\u2082 = '+T2.toFixed(2)+' N\nBelt Force (T\u2081+T\u2082) = '+(T1+T2).toFixed(2)+' N\nPulley Weight W = '+Wp+' N\nTotal Force = '+bf.toFixed(2)+' N\n\nR<sub>A</sub> = '+RAv.toFixed(2)+' N\nR<sub>B</sub> = '+RBv.toFixed(2)+' N'};
}

// ========== TWO GEARS (FIXED R_BV) ==========
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

    // VERTICAL: Taking moments about A:
    // R_BV x AB = F_tD(AD) - F_rC(AC)
    var V_RB=(FtD*AD-FrC*AC)/AB;
    var V_RA=FtD-FrC-V_RB;
    if(Wc>0||Wd>0){V_RB=((FtD+Wd)*AD-(FrC-Wc)*AC)/AB;V_RA=(FtD+Wd)-(FrC-Wc)-V_RB;}

    // HORIZONTAL: Taking moments about A:
    // R_BH x AB + F_tD(AD) = F_tC(AC)
    var H_RB=(FtC*AC-FtD*AD)/AB;
    var H_RB_neg=H_RB<0;
    var H_RA=FtC+H_RB-FtD;
    var H_RA_neg=H_RA<0;

    var M_CV=V_RA*AC;
    var M_DV=V_RB*(AB-AD);
    var M_CH=H_RA*AC;
    var M_DH=H_RB_neg?(-Math.abs(H_RB)*(AB-AD)):(H_RA*AD-FtC*(AD-AC));
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);
    var M_D=Math.sqrt(M_DV*M_DV+M_DH*M_DH);

    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);
        var vFC=(Wc>0)?(FrC-Wc):FrC;
        var vFD=(Wd>0)?(FtD+Wd):FtD;
        if(xi<AC){SFv.push(V_RA);BMv.push(V_RA*xi);SFh.push(H_RA);BMh.push(H_RA*xi);}
        else if(xi<AD){SFv.push(V_RA-vFC);BMv.push(V_RA*xi-vFC*(xi-AC));SFh.push(H_RA-FtC);BMh.push(H_RA*xi-FtC*(xi-AC));}
        else{SFv.push(V_RA-vFC+vFD);BMv.push(V_RA*xi-vFC*(xi-AC)+vFD*(xi-AD));SFh.push(H_RA-FtC+FtD);BMh.push(H_RA*xi-FtC*(xi-AC)+FtD*(xi-AD));}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));

    var det='';
    det+='Gear C: F<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = 2 \u00D7 '+T.toFixed(2)+' / '+Dc+' = '+FtC.toFixed(2)+' N\n';
    det+='F<sub>rC</sub> = F<sub>tC</sub> \u00D7 tan \u03B1 = '+FtC.toFixed(2)+' \u00D7 tan '+document.getElementById('twoGear_pressureAngle').value+'\u00B0 = '+FrC.toFixed(2)+' N\n\n';
    det+='Gear D: F<sub>tD</sub> = 2M<sub>t</sub>/D<sub>D</sub> = 2 \u00D7 '+T.toFixed(2)+' / '+Dd+' = '+FtD.toFixed(2)+' N\n';
    det+='F<sub>rD</sub> = F<sub>tD</sub> \u00D7 tan \u03B1 = '+FtD.toFixed(2)+' \u00D7 tan '+document.getElementById('twoGear_pressureAngle').value+'\u00B0 = '+FrD.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\n';
    det+='Taking moments about A:\nR<sub>BV</sub> \u00D7 '+AB+' = F<sub>tD</sub>(AD) \u2212 F<sub>rC</sub>(AC)\n';
    det+='R<sub>BV</sub> \u00D7 '+AB+' = ('+FtD.toFixed(2)+' \u00D7 '+AD+') \u2212 ('+FrC.toFixed(2)+' \u00D7 '+AC+')\n';
    det+='\u2234 R<sub>BV</sub> = '+V_RB.toFixed(2)+' N\n\n';
    det+='Also R<sub>AV</sub> + R<sub>BV</sub> = F<sub>tD</sub> \u2212 F<sub>rC</sub>\n';
    det+='\u2234 R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';
    det+='Moments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\n';
    det+='M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+M_DV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\n';
    det+='Taking moments about A:\nR<sub>BH</sub> \u00D7 '+AB+' + F<sub>tD</sub>(AD) = F<sub>tC</sub>(AC)\n';
    det+='\u2234 R<sub>BH</sub> = '+H_RB.toFixed(2)+' N';
    if(H_RB_neg)det+=' (acts opposite direction)';
    det+='\n\u2234 R<sub>AH</sub> = '+H_RA.toFixed(2)+' N';
    if(H_RA_neg)det+=' (acts opposite direction)';
    det+='\n\nMoments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DH</sub> = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\n';
    det+='M<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(0)+'\u00B2 + '+M_CH.toFixed(0)+'\u00B2) = '+M_C.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>D</sub> = \u221A(M\u00B2<sub>DV</sub> + M\u00B2<sub>DH</sub>) = \u221A('+M_DV.toFixed(0)+'\u00B2 + '+M_DH.toFixed(0)+'\u00B2) = '+M_D.toFixed(2)+' N\u00B7mm\n\n';
    det+='Max BM at '+(M_C>M_D?'C':'D')+': M = '+Math.max(M_C,M_D).toFixed(2)+' N\u00B7mm';

    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:Math.max.apply(null,BR),details:det};
}

// ========== PULLEY + GEAR ==========
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
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);
        if(xi<bC){SFv.push(-vP);BMv.push(-vP*xi);SFh.push(0);BMh.push(0);}
        else if(xi<bD){SFv.push(-vP+RCv);BMv.push(-vP*xi+RCv*(xi-bC));SFh.push(RCh);BMh.push(RCh*(xi-bC));}
        else if(xi<=gp){SFv.push(-vP+RCv+RDv);BMv.push(-vP*xi+RCv*(xi-bC)+RDv*(xi-bD));SFh.push(RCh+RDh);BMh.push(RCh*(xi-bC)+RDh*(xi-bD));}
        else{SFv.push(0);BMv.push(0);SFh.push(0);BMh.push(0);}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';
    det+='F<sub>tB</sub> = '+FtB.toFixed(2)+' N, F<sub>rB</sub> = '+FrB.toFixed(2)+' N\n\n';
    det+='Vertical:\nR<sub>CV</sub> = '+RCv.toFixed(2)+' N\nR<sub>DV</sub> = '+RDv.toFixed(2)+' N\n\n';
    det+='Horizontal:\nR<sub>CH</sub> = '+RCh.toFixed(2)+' N\nR<sub>DH</sub> = '+RDh.toFixed(2)+' N';
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:Math.max.apply(null,BR),details:det};
}

// ========== MULTI PULLEY ==========
function calcMultiPulley(T){
    var AB=parseFloat(document.getElementById('mp_bearingDist').value);
    var pA=parseFloat(document.getElementById('mp_pulleyA_pos').value);
    var WpA=parseFloat(document.getElementById('mp_pulleyA_weight').value);
    var T1A=parseFloat(document.getElementById('mp_T1A').value);var T2A=parseFloat(document.getElementById('mp_T2A').value);
    var pB=parseFloat(document.getElementById('mp_pulleyB_pos').value);
    var WpB=parseFloat(document.getElementById('mp_pulleyB_weight').value);
    var T1B=parseFloat(document.getElementById('mp_T1B').value);var T2B=parseFloat(document.getElementById('mp_T2B').value);
    var fA=T1A+T2A+WpA;var fB=T1B+T2B+WpB;
    var R2v=(fA*pA+fB*pB)/AB;var R1v=fA+fB-R2v;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);SFh.push(0);BMh.push(0);
        if(xi<pA){SFv.push(R1v);BMv.push(R1v*xi);}
        else if(xi<pB){SFv.push(R1v-fA);BMv.push(R1v*xi-fA*(xi-pA));}
        else{SFv.push(R1v-fA-fB);BMv.push(R1v*xi-fA*(xi-pA)-fB*(xi-pB));}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.abs(BMv[i]));
    var det='Pulley A: T\u2081='+T1A+'N, T\u2082='+T2A+'N, W='+WpA+'N\nTotal = '+fA.toFixed(2)+' N\n\n';
    det+='Pulley B: T\u2081='+T1B+'N, T\u2082='+T2B+'N, W='+WpB+'N\nTotal = '+fB.toFixed(2)+' N\n\n';
    det+='R\u2081 = '+R1v.toFixed(2)+' N\nR\u2082 = '+R2v.toFixed(2)+' N';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:Math.max.apply(null,BR),details:det};
}

// ========== DISPLAY RESULTS ==========
function displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText){
    document.getElementById('results').style.display='block';
    document.getElementById('diagrams').style.display='block';
    document.getElementById('torqueResult').textContent=T.toFixed(2)+' N\u00B7mm ('+(T/1000).toFixed(2)+' N\u00B7m)';
    document.getElementById('momentResult').textContent=maxBM.toFixed(2)+' N\u00B7mm';
    if(shaftType==='solid')document.getElementById('diameterResult').textContent='d = '+outerDia+' mm';
    else document.getElementById('diameterResult').textContent='d\u2092 = '+outerDia+' mm, d\u1D62 = '+innerDia.toFixed(2)+' mm';

    var h='<h3>Step-by-Step Solution</h3>';
    h+='<p><strong>Torque:</strong> T = (9.55 \u00D7 10\u2076 \u00D7 '+document.getElementById('power').value+') / '+document.getElementById('speed').value+' = '+T.toFixed(2)+' N\u00B7mm</p>';
    h+='<p><strong>\u03C3<sub>max</sub> = '+sigmaMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>\u03C4<sub>max</sub> = '+tauMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>C<sub>m</sub> = '+Cm.toFixed(2)+', C<sub>t</sub> = '+Ct.toFixed(2)+'</strong></p>';
    if(document.getElementById('hasKeyway').checked)h+='<p><strong>\u26A0 Keyway: stresses reduced by 25%</strong></p>';
    h+='<h3>Force Analysis</h3>';
    h+='<div style="background:#f0f0f0;padding:15px;border-radius:5px;line-height:1.8;font-size:14px;">'+result.details.replace(/\n/g,'<br>')+'</div>';
    h+='<h3>Maximum Bending Moment</h3>';
    h+='<p><strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong></p>';
    h+='<h3>Diameter Calculation</h3>';
    h+='<p><strong>Eq.(i) Max Normal Stress Theory ...3.6(a)/Pg 51, DHB:</strong></p>';
    h+='<p>d = [16/(\u03C0\u03C3<sub>max</sub>) \u00D7 {C<sub>m</sub>M + \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)}]<sup>1/3</sup></p>';
    h+='<p>= [16/(\u03C0 \u00D7 '+sigmaMax.toFixed(2)+') \u00D7 {'+Cm+' \u00D7 '+maxBM.toFixed(0)+' + \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(0)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(0)+')\u00B2)}]<sup>1/3</sup></p>';
    h+='<p><strong>\u2234 d = '+dNormal.toFixed(2)+' mm ...Eq.(i)</strong></p>';
    h+='<br><p><strong>Eq.(ii) Max Shear Stress Theory ...3.6(b)/Pg 51, DHB:</strong></p>';
    h+='<p>d = [16/(\u03C0\u03C4<sub>max</sub>) \u00D7 \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)]<sup>1/3</sup></p>';
    h+='<p>= [16/(\u03C0 \u00D7 '+tauMax.toFixed(2)+') \u00D7 \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(0)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(0)+')\u00B2)]<sup>1/3</sup></p>';
    h+='<p><strong>\u2234 d = '+dShear.toFixed(2)+' mm ...Eq.(ii)</strong></p>';
    h+='<br><p>Based on Eqs. (i) and (ii), select maximum diameter:</p>';
    h+='<p>d = '+Math.max(dNormal,dShear).toFixed(2)+' mm</p>';
    h+='<p style="color:red;font-size:1.3em;font-weight:bold;">\u2234 Standard size of shaft, d = '+outerDia+' mm ...Tb. 3.5(a)/Pg 57, DHB</p>';
    if(deflText){h+='<h3>Angular Deflection \u03B8</h3>'+deflText;}
    document.getElementById('detailedResults').innerHTML=h;
    document.getElementById('results').scrollIntoView({behavior:'smooth'});
}

// ========== GENERATE
