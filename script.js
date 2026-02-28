function toggleProblemType(){
    var v=document.querySelector('input[name="problemType"]:checked').value;
    var ids=['normalSection','udlSection','singlePulleySection','twoGearSection','pulleyGearSection','twoPulleyVHSection','gearAngledPulleySection'];
    var map={normal:'normalSection',udl:'udlSection',singlepulley:'singlePulleySection',twogear:'twoGearSection',pulleygear:'pulleyGearSection',twopulleyvh:'twoPulleyVHSection',gearangledpulley:'gearAngledPulleySection'};
    for(var i=0;i<ids.length;i++){document.getElementById(ids[i]).style.display='none';}
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
    document.getElementById('frictionTensionInputs').style.display=v==='friction'?'block':'none';
}
function toggleSPBearings(){
    var v=document.querySelector('input[name="sp_bearingCount"]:checked').value;
    document.getElementById('sp_1bearing').style.display=v==='1'?'block':'none';
    document.getElementById('sp_2bearings').style.display=v==='2'?'block':'none';
}
function toggleGearDiaInput(){
    var v=document.querySelector('input[name="gearDiaInput"]:checked').value;
    document.getElementById('gearC_directDia').style.display=v==='direct'?'block':'none';
    document.getElementById('gearC_moduleDia').style.display=v==='module'?'block':'none';
    document.getElementById('gearD_directDia').style.display=v==='direct'?'block':'none';
    document.getElementById('gearD_moduleDia').style.display=v==='module'?'block':'none';
}
function toggleVHC(){var v=document.querySelector('input[name="vhC_tensionMethod"]:checked').value;document.getElementById('vhC_directInputs').style.display=v==='direct'?'block':'none';document.getElementById('vhC_frictionInputs').style.display=v==='friction'?'block':'none';}
function toggleVHD(){var v=document.querySelector('input[name="vhD_tensionMethod"]:checked').value;document.getElementById('vhD_directInputs').style.display=v==='direct'?'block':'none';document.getElementById('vhD_frictionInputs').style.display=v==='friction'?'block':'none';}
document.addEventListener('DOMContentLoaded',function(){document.getElementById('calcDeflection').addEventListener('change',function(){document.getElementById('deflectionInputs').style.display=this.checked?'block':'none';});});
function evalFraction(s){s=s.trim();if(s.indexOf('/')!==-1){var p=s.split('/');return parseFloat(p[0])/parseFloat(p[1]);}return parseFloat(s);}

function calculateShaft(){
    var problemType=document.querySelector('input[name="problemType"]:checked').value;
    var P=parseFloat(document.getElementById('power').value);
    var N=parseFloat(document.getElementById('speed').value);
    var T=(9.55e6*P)/N;
    var stressMethod=document.querySelector('input[name="stressMethod"]:checked').value;
    var sigmaMax,tauMax;
    if(stressMethod==='direct'){tauMax=parseFloat(document.getElementById('tauAllow').value);sigmaMax=2*tauMax;}
    else if(stressMethod==='material'){var sigUt=parseFloat(document.getElementById('sigmaUt').value);var sigYt=parseFloat(document.getElementById('sigmaYt').value);var nP=parseFloat(document.getElementById('nPrime').value);var sigU=sigUt/nP;var sigY=sigYt/nP;sigmaMax=Math.min(0.36*sigU,0.6*sigY);tauMax=Math.min(0.18*sigU,0.3*sigY);}
    else if(stressMethod==='working'){sigmaMax=parseFloat(document.getElementById('sigmaWorking').value);tauMax=parseFloat(document.getElementById('tauWorking').value);}
    else if(stressMethod==='yield'){var sYt=parseFloat(document.getElementById('sigmaYtYield').value);var fos=parseFloat(document.getElementById('fosYield').value);sigmaMax=sYt/fos;tauMax=sigmaMax/2;}
    if(document.getElementById('hasKeyway').checked){sigmaMax=0.75*sigmaMax;tauMax=0.75*tauMax;}
    var Cm=parseFloat(document.getElementById('cm').value);var Ct=parseFloat(document.getElementById('ct').value);
    var result;
    if(problemType==='normal') result=calcNormalShaft(T);
    else if(problemType==='udl') result=calcUDL(T);
    else if(problemType==='singlepulley') result=calcSinglePulley(T);
    else if(problemType==='twogear') result=calcTwoGear(T);
    else if(problemType==='pulleygear') result=calcPulleyGear(T);
    else if(problemType==='twopulleyvh') result=calcTwoPulleyVH(T);
    else if(problemType==='gearangledpulley') result=calcGearAngledPulley(T);
    if((problemType==='twopulleyvh'||problemType==='pulleygear')&&result.torqueOverride!==undefined){T=result.torqueOverride;}
    var shaftType=document.querySelector('input[name="shaftType"]:checked').value;
    var maxBM=result.maxBM;
    var useShearOnly=(problemType==='twopulleyvh'||problemType==='pulleygear'||problemType==='singlepulley');
    var dNormal,dShear;
    if(useShearOnly){var t2s=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));dShear=Math.pow((16/(Math.PI*tauMax))*t2s,1/3);dNormal=dShear;}
    else{var t1=Cm*maxBM;var t2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));dNormal=Math.pow((16/(Math.PI*sigmaMax))*(t1+t2),1/3);dShear=Math.pow((16/(Math.PI*tauMax))*t2,1/3);}
    var outerDia,innerDia;
    var stdSizes=[6,8,10,12,14,16,18,20,22,25,28,32,36,40,45,50,56,63,71,80,90,100,110,125,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,450,480,500,530,560,600];
    if(shaftType==='solid'){var dd=Math.max(dNormal,dShear);outerDia=null;for(var i=0;i<stdSizes.length;i++){if(stdSizes[i]>=dd){outerDia=stdSizes[i];break;}}if(!outerDia)outerDia=Math.ceil(dd/10)*10;innerDia=0;}
    else{var k=evalFraction(document.getElementById('diameterRatio').value);if(useShearOnly){var t2h=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));var dd2=Math.pow((16/(Math.PI*tauMax))*t2h*(1/(1-Math.pow(k,4))),1/3);}else{var t1h=Cm*maxBM;var t2h2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));var doN=Math.pow((16/(Math.PI*sigmaMax))*(t1h+t2h2)*(1/(1-Math.pow(k,4))),1/3);var doS=Math.pow((16/(Math.PI*tauMax))*t2h2*(1/(1-Math.pow(k,4))),1/3);var dd2=Math.max(doN,doS);}outerDia=null;for(var i=0;i<stdSizes.length;i++){if(stdSizes[i]>=dd2){outerDia=stdSizes[i];break;}}if(!outerDia)outerDia=Math.ceil(dd2/10)*10;innerDia=k*outerDia;}
    var deflText='';
    if(document.getElementById('calcDeflection').checked){var Gd=parseFloat(document.getElementById('deflG').value)*1000;var Ld=parseFloat(document.getElementById('deflL').value);var D4=Math.pow(outerDia,4);var theta=(584*T*Ld)/(Gd*D4);deflText='<p>\u03B8 = (584 \u00D7 T \u00D7 L) / (G \u00D7 d\u2074)</p><p>= (584 \u00D7 '+T.toFixed(2)+' \u00D7 '+Ld+') / ('+Gd+' \u00D7 '+outerDia+'\u2074)</p><p><strong>\u03B8 = '+theta.toFixed(4)+' degrees</strong></p>';}
    displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText,problemType);
    generateDiagrams(result,T);
}

function calcNormalShaft(T){
    var AB=parseFloat(document.getElementById('normal_bearingDist').value);var gp=parseFloat(document.getElementById('normal_gearPos').value);var gd=parseFloat(document.getElementById('normal_gearDia').value);var alpha=parseFloat(document.getElementById('normal_pressureAngle').value)*Math.PI/180;var Wg=parseFloat(document.getElementById('normal_gearWeight').value);
    var Ft=(2*T)/gd;var Fr=Ft*Math.tan(alpha);var RAh=Ft*(AB-gp)/AB;var RBh=Ft-RAh;var vf=Fr+Wg;var RAv=vf*(AB-gp)/AB;var RBv=vf-RAv;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);if(xi<gp){SFh.push(RAh);BMh.push(RAh*xi);SFv.push(RAv);BMv.push(RAv*xi);}else{SFh.push(RAh-Ft);BMh.push(RAh*xi-Ft*(xi-gp));SFv.push(RAv-vf);BMv.push(RAv*xi-vf*(xi-gp));}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));var maxBM=0;for(var i=0;i<n;i++){if(BR[i]>maxBM)maxBM=BR[i];}
    var det='F<sub>t</sub> = 2T/D = '+Ft.toFixed(2)+' N\nF<sub>r</sub> = F<sub>t</sub> \u00D7 tan \u03B1 = '+Fr.toFixed(2)+' N\n\nR<sub>AH</sub> = '+RAh.toFixed(2)+' N, R<sub>BH</sub> = '+RBh.toFixed(2)+' N\nR<sub>AV</sub> = '+RAv.toFixed(2)+' N, R<sub>BV</sub> = '+RBv.toFixed(2)+' N';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcUDL(T){
    var L=parseFloat(document.getElementById('totalLength').value);var lU=parseFloat(document.getElementById('udlLength').value);var w=parseFloat(document.getElementById('udlIntensity').value);var b1=(L-lU)/2;var RA=(w*lU)/2;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);if(xi<b1||xi>b1+lU){SFh.push(0);BMh.push(0);SFv.push(0);BMv.push(0);}else{var d=xi-b1;SFh.push(RA-w*d);BMh.push(RA*d-w*d*d/2);SFv.push(RA-w*d);BMv.push(RA*d-w*d*d/2);}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));var maxBM=0;for(var i=0;i<n;i++){if(BR[i]>maxBM)maxBM=BR[i];}
    var det='w = '+w+' N/mm, l = '+lU+' mm\nR<sub>A</sub> = R<sub>B</sub> = '+RA.toFixed(2)+' N\nM<sub>max</sub> = '+(w*lU*lU/8).toFixed(2)+' N\u00B7mm';
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcSinglePulley(T){var bc=document.querySelector('input[name="sp_bearingCount"]:checked').value;return(bc==='1')?calcSinglePulley1B(T):calcSinglePulley2B(T);}

function calcSinglePulley1B(T){
    var L=parseFloat(document.getElementById('sp_overhangDist').value);var pd=parseFloat(document.getElementById('sp_pulleyDia1').value);var Wp=parseFloat(document.getElementById('sp_pulleyWeight1').value);var R=pd/2;
    var T1,T2,td='';var tM=document.querySelector('input[name="tensionInput"]:checked').value;
    if(tM==='direct'){T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);td='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';}
    else if(tM==='ratio'){var r=parseFloat(document.getElementById('sp_tensionRatio').value);T2=T/(R*(r-1));T1=r*T2;td='T\u2081/T\u2082 = '+r+'\nT = (T\u2081 \u2212 T\u2082) \u00D7 R = '+T.toFixed(2)+'\n(T\u2081 \u2212 T\u2082) = '+(T/R).toFixed(2)+' N\nT\u2082 = '+T2.toFixed(2)+' N, T\u2081 = '+T1.toFixed(2)+' N\n';}
    else{var mu=parseFloat(document.getElementById('sp_mu').value);var tr=parseFloat(document.getElementById('sp_theta').value)*Math.PI/180;var ratio=Math.exp(mu*tr);T2=T/(R*(ratio-1));T1=ratio*T2;td='T\u2081/T\u2082 = e^(\u03BC\u03B8) = '+ratio.toFixed(4)+'\nT\u2082 = '+T2.toFixed(2)+' N, T\u2081 = '+T1.toFixed(2)+' N\n';}
    var tf=T1+T2+Wp;var maxBM=tf*L;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);SFh.push(0);BMh.push(0);SFv.push(tf);BMv.push(tf*xi);}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.abs(BMv[i]));
    var det='\u2550\u2550\u2550 1 Bearing (Cantilever) \u2550\u2550\u2550\n\n'+td+'\nM = (T\u2081 + T\u2082 + W) \u00D7 L = '+tf.toFixed(2)+' \u00D7 '+L+'\n<strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcSinglePulley2B(T){
    var AB=parseFloat(document.getElementById('sp_bearingDist').value);var AC=parseFloat(document.getElementById('sp_pulleyPos').value);var pd=parseFloat(document.getElementById('sp_pulleyDia2').value);var Wp=parseFloat(document.getElementById('sp_pulleyWeight2').value);var R=pd/2;
    var T1,T2,td='';var tM=document.querySelector('input[name="tensionInput"]:checked').value;
    if(tM==='direct'){T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);td='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';}
    else if(tM==='ratio'){var r=parseFloat(document.getElementById('sp_tensionRatio').value);T2=T/(R*(r-1));T1=r*T2;td='T = (T\u2081 \u2212 T\u2082) \u00D7 R\nT\u2081 \u2212 T\u2082 = '+(T/R).toFixed(2)+' N\nT\u2081 = '+r+' \u00D7 T\u2082\nT\u2082 = '+T2.toFixed(2)+' N, T\u2081 = '+T1.toFixed(2)+' N\n';}
    else{var mu=parseFloat(document.getElementById('sp_mu').value);var tr=parseFloat(document.getElementById('sp_theta').value)*Math.PI/180;var ratio=Math.exp(mu*tr);T2=T/(R*(ratio-1));T1=ratio*T2;td='T\u2081/T\u2082 = e^(\u03BC\u03B8) = '+ratio.toFixed(4)+'\nT\u2082 = '+T2.toFixed(2)+' N, T\u2081 = '+T1.toFixed(2)+' N\n';}
    var RBV=(Wp*AC)/AB;var RAV=Wp-RBV;var Fh=T1+T2;var RBH=(Fh*AC)/AB;var RAH=Fh-RBH;
    var MCV=RAV*AC;var MCH=RAH*AC;var MC=Math.sqrt(MCV*MCV+MCH*MCH);
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);if(xi<AC){SFv.push(RAV);BMv.push(RAV*xi);SFh.push(RAH);BMh.push(RAH*xi);}else{SFv.push(RAV-Wp);BMv.push(RAV*xi-Wp*(xi-AC));SFh.push(RAH-Fh);BMh.push(RAH*xi-Fh*(xi-AC));}}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 2 Bearings \u2550\u2550\u2550\n\n'+td+'\n\u2550\u2550\u2550 Vertical \u2550\u2550\u2550\nR<sub>BV</sub> = '+RBV.toFixed(2)+' N, R<sub>AV</sub> = '+RAV.toFixed(2)+' N\nM<sub>CV</sub> = '+MCV.toFixed(2)+' N\u00B7mm\n\n\u2550\u2550\u2550 Horizontal \u2550\u2550\u2550\nR<sub>BH</sub> = '+RBH.toFixed(2)+' N, R<sub>AH</sub> = '+RAH.toFixed(2)+' N\nM<sub>CH</sub> = '+MCH.toFixed(2)+' N\u00B7mm\n\n\u2550\u2550\u2550 Resultant \u2550\u2550\u2550\n<strong>M<sub>C</sub> = '+MC.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:MC,details:det};
}

function calcTwoGear(T){
    var AB=parseFloat(document.getElementById('bearingDistance').value);var AC=parseFloat(document.getElementById('gearC_pos').value);var AD=parseFloat(document.getElementById('gearD_pos').value);
    var diaInput=document.querySelector('input[name="gearDiaInput"]:checked').value;var Dc,Dd,diaDetail='';
    if(diaInput==='direct'){Dc=parseFloat(document.getElementById('gearC_dia').value);Dd=parseFloat(document.getElementById('gearD_dia').value);}
    else{var mC=parseFloat(document.getElementById('gearC_module').value);var zC=parseFloat(document.getElementById('gearC_teeth').value);var mD=parseFloat(document.getElementById('gearD_module').value);var zD=parseFloat(document.getElementById('gearD_teeth').value);Dc=mC*zC;Dd=mD*zD;diaDetail='D<sub>C</sub> = m \u00D7 z = '+mC+' \u00D7 '+zC+' = '+Dc+' mm\nD<sub>D</sub> = m \u00D7 z = '+mD+' \u00D7 '+zD+' = '+Dd+' mm\n\n';}
    var Wc=parseFloat(document.getElementById('gearC_weight').value);var Wd=parseFloat(document.getElementById('gearD_weight').value);var alpha=parseFloat(document.getElementById('twoGear_pressureAngle').value)*Math.PI/180;
    var fc=document.querySelector('input[name="twoGearForceConfig"]:checked').value;
    var FtC=(2*T)/Dc;var FrC=FtC*Math.tan(alpha);var FtD=(2*T)/Dd;var FrD=FtD*Math.tan(alpha);
    var V_RB,V_RA,H_RB,H_RA,M_CV,M_DV,M_CH,M_DH;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    var det='';if(diaDetail)det+=diaDetail;
    det+='F<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = '+FtC.toFixed(2)+' N\nF<sub>rC</sub> = '+FrC.toFixed(2)+' N\nF<sub>tD</sub> = 2M<sub>t</sub>/D<sub>D</sub> = '+FtD.toFixed(2)+' N\nF<sub>rD</sub> = '+FrD.toFixed(2)+' N\n\n';

    if(fc==='crossed'){
        V_RB=(FtD*AD-FrC*AC)/AB;V_RA=FtD-FrC-V_RB;
        var HRBr=(FtC*AC-FrD*AD)/AB;var HRBn=HRBr<0;var HRBa=Math.abs(HRBr);H_RB=HRBr;
        if(HRBn){H_RA=FtC+HRBa-FrD;}else{H_RA=FtC-FrD-H_RB;}
        M_CV=V_RA*AC;M_DV=V_RB*(AB-AD);M_CH=H_RA*AC;M_DH=H_RB*(AB-AD);
        for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);
            if(xi<AC){SFv.push(V_RA);BMv.push(V_RA*xi);SFh.push(H_RA);BMh.push(H_RA*xi);}
            else if(xi<AD){SFv.push(V_RA+FrC);BMv.push(V_RA*xi+FrC*(xi-AC));SFh.push(H_RA-FtC);BMh.push(H_RA*xi-FtC*(xi-AC));}
            else{SFv.push(V_RA+FrC-FtD);BMv.push(V_RA*xi+FrC*(xi-AC)-FtD*(xi-AD));SFh.push(H_RA-FtC+FrD);BMh.push(H_RA*xi-FtC*(xi-AC)+FrD*(xi-AD));}}
        det+='\u2550\u2550\u2550 Vertical (Crossed) \u2550\u2550\u2550\nR<sub>BV</sub>(AB) = F<sub>tD</sub>(AD) \u2212 F<sub>rC</sub>(AC)\nR<sub>BV</sub> = '+V_RB.toFixed(2)+' N\nR<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';
        det+='\u2550\u2550\u2550 Horizontal (Crossed) \u2550\u2550\u2550\nR<sub>BH</sub>(AB) + F<sub>rD</sub>(AD) = F<sub>tC</sub>(AC)\nR<sub>BH</sub> = '+HRBr.toFixed(2)+' N';
        if(HRBn)det+=' (opposite direction, R<sub>BH</sub> = '+HRBa.toFixed(2)+' N)';
        det+='\nR<sub>AH</sub> = '+H_RA.toFixed(2)+' N\n\n';
    }else{
        var vC=FtC+Wc;var vD=FtD+Wd;V_RB=(vC*AC+vD*AD)/AB;V_RA=vC+vD-V_RB;
        var HRBr=(FrC*AC-FrD*AD)/AB;var HRBn=HRBr<0;var HRBa=Math.abs(HRBr);H_RB=HRBr;
        if(HRBn){H_RA=FrC+HRBa-FrD;}else{H_RA=FrC+H_RB-FrD;}
        M_CV=V_RA*AC;M_DV=V_RB*(AB-AD);M_CH=H_RA*AC;M_DH=H_RB*(AB-AD);
        var hA=(H_RA<0)?-Math.abs(H_RA):H_RA;
        for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);
            if(xi<AC){SFv.push(V_RA);BMv.push(V_RA*xi);SFh.push(hA);BMh.push(hA*xi);}
            else if(xi<AD){SFv.push(V_RA-vC);BMv.push(V_RA*xi-vC*(xi-AC));SFh.push(hA-FrC);BMh.push(hA*xi-FrC*(xi-AC));}
            else{SFv.push(V_RA-vC+vD);BMv.push(V_RA*xi-vC*(xi-AC)+vD*(xi-AD));SFh.push(hA-FrC+FrD);BMh.push(hA*xi-FrC*(xi-AC)+FrD*(xi-AD));}}
        det+='\u2550\u2550\u2550 Vertical (Same) \u2550\u2550\u2550\nR<sub>BV</sub> = '+V_RB.toFixed(2)+' N, R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';
        det+='\u2550\u2550\u2550 Horizontal (Same) \u2550\u2550\u2550\nR<sub>BH</sub> = '+HRBr.toFixed(2)+' N';if(HRBn)det+=' (opposite)';det+='\nR<sub>AH</sub> = '+H_RA.toFixed(2)+' N\n\n';
    }
    det+='Moments:\nM<sub>CV</sub> = '+M_CV.toFixed(2)+' N\u00B7mm, M<sub>DV</sub> = '+M_DV.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>CH</sub> = '+M_CH.toFixed(2)+' N\u00B7mm, M<sub>DH</sub> = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';
    var MC=Math.sqrt(M_CV*M_CV+M_CH*M_CH);var MD=Math.sqrt(M_DV*M_DV+M_DH*M_DH);var maxBM=Math.max(MC,MD);
    det+='\u2550\u2550\u2550 Resultant \u2550\u2550\u2550\nM<sub>C</sub> = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2) = '+MC.toFixed(2)+' N\u00B7mm\nM<sub>D</sub> = \u221A('+M_DV.toFixed(2)+'\u00B2 + ('+M_DH.toFixed(2)+')\u00B2) = '+MD.toFixed(2)+' N\u00B7mm\n\n';
    det+='<strong>Max BM at '+(MC>MD?'C':'D')+': M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    var BR=[];for(var i=0;i<x.length;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcPulleyGear(T){
    var AB=parseFloat(document.getElementById('pg_bearingDistance').value);var AC=parseFloat(document.getElementById('pulley_pos').value);var pd=parseFloat(document.getElementById('pulley_dia').value);var Wa=parseFloat(document.getElementById('pulley_weight').value);var tr=parseFloat(document.getElementById('tensionRatio').value);var AD=parseFloat(document.getElementById('gear_pos').value);var gd=parseFloat(document.getElementById('gear_dia').value);var Wb=parseFloat(document.getElementById('gear_weight').value);var alpha=parseFloat(document.getElementById('pg_pressureAngle').value)*Math.PI/180;
    var Rp=pd/2;var T2v=T/(Rp*(tr-1));var T1v=tr*T2v;var FtD=(2*T)/gd;var FrD=FtD*Math.tan(alpha);
    var FvD=FtD+Wb;var RBV=(FvD*AD)/AB;var RAV=FvD-RBV;var FhC=T1v+T2v+Wa;var RBH=(FrD*AD-FhC*AC)/AB;var RAH=FrD-FhC-RBH;
    var MCV=RAV*AC;var MDV=RBV*(AB-AD);var MCH=RAH*AC;var MDH=RBH*(AB-AD);
    var MC=Math.sqrt(MCV*MCV+MCH*MCH);var MD=Math.sqrt(MDV*MDV+MDH*MDH);var maxBM=Math.max(MC,MD);var mL=(MC>=MD)?'C':'D';
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);var sv=RAV,mv=RAV*xi,sh=RAH,mh=RAH*xi;if(xi>=AD){sv-=FvD;mv-=FvD*(xi-AD);}if(xi>=AC){sh-=FhC;mh-=FhC*(xi-AC);}if(xi>=AD){sh+=FrD;mh+=FrD*(xi-AD);}SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='T\u2081 = '+T1v.toFixed(2)+' N, T\u2082 = '+T2v.toFixed(2)+' N\nF<sub>tD</sub> = '+FtD.toFixed(2)+' N, F<sub>rD</sub> = '+FrD.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Vertical \u2550\u2550\u2550\nR<sub>BV</sub> = '+RBV.toFixed(2)+' N, R<sub>AV</sub> = '+RAV.toFixed(2)+' N\nM<sub>CV</sub> = '+MCV.toFixed(2)+', M<sub>DV</sub> = '+MDV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal \u2550\u2550\u2550\nR<sub>BH</sub> = '+RBH.toFixed(2)+' N, R<sub>AH</sub> = '+RAH.toFixed(2)+' N';if(RAH<0)det+=' (opposite)';
    det+='\nM<sub>CH</sub> = '+MCH.toFixed(2)+', M<sub>DH</sub> = '+MDH.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Resultant \u2550\u2550\u2550\nM<sub>C</sub> = '+MC.toFixed(2)+' N\u00B7mm\nM<sub>D</sub> = '+MD.toFixed(2)+' N\u00B7mm\n\n<strong>Max BM at '+mL+': M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det,torqueOverride:T};
}

function calcTwoPulleyVH(Ta){
    var AB=parseFloat(document.getElementById('vh_bearingDist').value);var pC=parseFloat(document.getElementById('vh_pulleyC_pos').value);var dC=parseFloat(document.getElementById('vh_pulleyC_dia').value);var pD=parseFloat(document.getElementById('vh_pulleyD_pos').value);var dD=parseFloat(document.getElementById('vh_pulleyD_dia').value);var RC=dC/2;var RD=dD/2;
    var T3,T4,dMS='';var dM=document.querySelector('input[name="vhD_tensionMethod"]:checked').value;
    if(dM==='direct'){T3=parseFloat(document.getElementById('vh_T3D').value);T4=parseFloat(document.getElementById('vh_T4D').value);}else{var T3f=parseFloat(document.getElementById('vh_T3D_f').value);var muD=parseFloat(document.getElementById('vh_muD').value);var tDr=parseFloat(document.getElementById('vh_thetaD').value)*Math.PI/180;var rD=Math.exp(muD*tDr);T3=T3f;T4=T3/rD;dMS='T\u2083/T\u2084 = '+rD.toFixed(4)+', T\u2084 = '+T4.toFixed(2)+' N\n';}
    var TR=(T3-T4)*RD;var T1,T2,cMS='';var cM=document.querySelector('input[name="vhC_tensionMethod"]:checked').value;
    if(cM==='direct'){T1=parseFloat(document.getElementById('vh_T1C').value);T2=parseFloat(document.getElementById('vh_T2C').value);}else{var muC=parseFloat(document.getElementById('vh_muC').value);var tCr=parseFloat(document.getElementById('vh_thetaC').value)*Math.PI/180;var rC=Math.exp(muC*tCr);T2=TR/(RC*(rC-1));T1=rC*T2;cMS='T\u2081/T\u2082 = '+rC.toFixed(4)+', T\u2082 = '+T2.toFixed(2)+' N, T\u2081 = '+T1.toFixed(2)+' N\n';}
    var tM2=document.querySelector('input[name="vhTorqueMethod"]:checked').value;var TU=(tM2==='pulley')?TR:Ta;
    var Fv=T1+T2;var RBV=(Fv*pC)/AB;var RAV=Fv-RBV;var Fh=T3+T4;var RBH=(Fh*pD)/AB;var RAH=Fh-RBH;
    var MCV=RAV*pC;var MCH=RAH*pC;var MBH=(pD>AB)?-Fh*(pD-AB):0;
    var MC=Math.sqrt(MCV*MCV+MCH*MCH);var MB=Math.abs(MBH);var maxBM=Math.max(MC,MB);var mL=(MC>=MB)?'C':'B';
    var sX=Math.min(0,pC)-10;var eX=Math.max(AB,pD)+10;var tL=eX-sX;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=sX+(i/(n-1))*tL;x.push(xi);var sv=0,mv=0,sh=0,mh=0;if(xi>=0){sv+=RAV;mv+=RAV*xi;}if(xi>=pC){sv-=Fv;mv-=Fv*(xi-pC);}if(xi>=AB){sv+=RBV;mv+=RBV*(xi-AB);}if(xi>=0){sh+=RAH;mh+=RAH*xi;}if(xi>=AB){sh+=RBH;mh+=RBH*(xi-AB);}if(xi>=pD){sh-=Fh;mh-=Fh*(xi-pD);}SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 Pulley D \u2550\u2550\u2550\n';if(dM==='friction')det+=dMS;else det+='T\u2083 = '+T3.toFixed(2)+', T\u2084 = '+T4.toFixed(2)+' N\n';
    det+='T<sub>R</sub> = '+TR.toFixed(2)+' N\u00B7mm\n\n\u2550\u2550\u2550 Pulley C \u2550\u2550\u2550\n';if(cM==='friction')det+=cMS;else det+='T\u2081 = '+T1.toFixed(2)+', T\u2082 = '+T2.toFixed(2)+' N\n';
    det+='\n\u2550\u2550\u2550 Vertical \u2550\u2550\u2550\nR<sub>BV</sub> = '+RBV.toFixed(2)+' N, R<sub>AV</sub> = '+RAV.toFixed(2)+' N\nM<sub>CV</sub> = '+MCV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal \u2550\u2550\u2550\nR<sub>BH</sub> = '+RBH.toFixed(2)+' N, R<sub>AH</sub> = '+RAH.toFixed(2)+' N';if(RAH<0)det+=' (opposite)';
    det+='\nM<sub>CH</sub> = '+MCH.toFixed(2)+' N\u00B7mm\n';if(pD>AB)det+='M<sub>BH</sub> = '+MBH.toFixed(2)+' N\u00B7mm\n';
    det+='\n\u2550\u2550\u2550 Resultant \u2550\u2550\u2550\nM<sub>C</sub> = '+MC.toFixed(2)+' N\u00B7mm\nM<sub>B</sub> = '+MB.toFixed(2)+' N\u00B7mm\n\n<strong>Max BM at '+mL+': M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:tL,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det,torqueOverride:TU};
}

// ========== GEAR + ANGLED PULLEY — FIXED to match textbook exactly ==========
// Layout: Bearings A (left) and B (right), distance AB
// Gear at C (between bearings): F_tC and W_C act vertical, F_rC acts horizontal
// Pulley at D (overhangs past B): belt at angle theta
//   Vertical at D: (T1+T2)*sin(theta) + W_D
//   Horizontal at D: (T1+T2)*cos(theta)
//
// VERTICAL (taking moments about A):
//   R_BV * AB = [F_tC + W_C] * AC + [(T1+T2)*sin(theta) + W_D] * AD
//   R_AV + R_BV = [F_tC + W_C] + [(T1+T2)*sin(theta) + W_D]
//   M_AV = M_DV = 0
//   M_CV = R_AV * AC
//   M_BV = -[(T1+T2)*sin(theta) + W_D] * (AD - AB)  (overhang moment at B)
//
// HORIZONTAL (taking moments about A):
//   R_BH * AB = F_rC * AC + [(T1+T2)*cos(theta)] * AD
//   R_AH + R_BH = F_rC + [(T1+T2)*cos(theta)]
//   M_AH = M_DH = 0
//   M_CH = R_AH * AC
//   M_BH = -[(T1+T2)*cos(theta)] * (AD - AB)  (overhang moment at B)
//
// Resultant: M_C = sqrt(M_CV^2 + M_CH^2), M_B = sqrt(M_BV^2 + M_BH^2)
function calcGearAngledPulley(T){
    var AB=parseFloat(document.getElementById('gap_bearingDist').value);
    var AC=parseFloat(document.getElementById('gap_gearPos').value);
    var gd=parseFloat(document.getElementById('gap_gearDia').value);
    var alphaDeg=document.getElementById('gap_pressureAngle').value;
    var alpha=parseFloat(alphaDeg)*Math.PI/180;
    var Wc=parseFloat(document.getElementById('gap_gearWeight').value);
    var AD=parseFloat(document.getElementById('gap_pulleyPos').value);
    var pDia=parseFloat(document.getElementById('gap_pulleyDia').value);
    var thetaDeg=parseFloat(document.getElementById('gap_beltAngle').value);
    var thetaRad=thetaDeg*Math.PI/180;
    var Wd=parseFloat(document.getElementById('gap_pulleyWeight').value);
    var tr=parseFloat(document.getElementById('gap_tensionRatio').value);

    var Rp=pDia/2;
    // T1/T2 = tr, T = (T1-T2)*Rp
    var T2val=T/(Rp*(tr-1));
    var T1val=tr*T2val;

    // Belt force components
    var beltSum=T1val+T2val;
    var beltV=beltSum*Math.sin(thetaRad);
    var beltH=beltSum*Math.cos(thetaRad);

    // Gear forces
    var FtC=(2*T)/gd;
    var FrC=FtC*Math.tan(alpha);

    // Total forces at C and D
    var Fv_C=FtC+Wc;       // vertical at C (downward)
    var Fv_D=beltV+Wd;     // vertical at D (downward)
    var Fh_D=beltH;        // horizontal at D

    // ===== VERTICAL PLANE =====
    // Taking moments about A:
    // R_BV * AB = Fv_C * AC + Fv_D * AD
    var R_BV=(Fv_C*AC+Fv_D*AD)/AB;
    // R_AV + R_BV = Fv_C + Fv_D
    var R_AV=Fv_C+Fv_D-R_BV;

    // Moments in vertical plane
    var M_AV=0;
    var M_DV=0;
    var M_CV=R_AV*AC;
    // M_BV: moment at B from overhang (if D overhangs past B)
    var M_BV=0;
    if(AD>AB){
        M_BV=-(Fv_D)*(AD-AB);
    }

    // ===== HORIZONTAL PLANE =====
    // Taking moments about A:
    // R_BH * AB = F_rC * AC + beltH * AD
    var R_BH=(FrC*AC+beltH*AD)/AB;
    // R_AH + R_BH = F_rC + beltH
    var R_AH=FrC+beltH-R_BH;
    var R_AH_neg=R_AH<0;

    // Moments in horizontal plane
    var M_AH=0;
    var M_DH=0;
    var M_CH=R_AH*AC;
    // M_BH: moment at B from overhang
    var M_BH=0;
    if(AD>AB){
        M_BH=-(beltH)*(AD-AB);
    }

    // Resultant moments
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);
    var M_B=Math.sqrt(M_BV*M_BV+M_BH*M_BH);
    var maxBM=Math.max(M_C,M_B);
    var maxLoc=(M_C>=M_B)?'C':'B';

    // SFD/BMD
    var startX=-10;
    var endX=Math.max(AB,AD)+10;
    var totalL=endX-startX;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=startX+(i/(n-1))*totalL;x.push(xi);
        var sv=0,mv=0,sh=0,mh=0;
        if(xi>=0){sv+=R_AV;mv+=R_AV*xi;}
        if(xi>=AC){sv-=Fv_C;mv-=Fv_C*(xi-AC);}
        if(xi>=AB){sv+=R_BV;mv+=R_BV*(xi-AB);}
        if(xi>=AD){sv-=Fv_D;mv-=Fv_D*(xi-AD);}
        if(xi>=0){sh+=R_AH;mh+=R_AH*xi;}
        if(xi>=AC){sh-=FrC;mh-=FrC*(xi-AC);}
        if(xi>=AB){sh+=R_BH;mh+=R_BH*(xi-AB);}
        if(xi>=AD){sh-=beltH;mh-=beltH*(xi-AD);}
        SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);
    }
    var BR=[];for(var i=0;i<n;i++)BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));

    // Build detailed solution matching textbook step by step
    var det='';
    det+='\u2550\u2550\u2550 Analysis of Pulley D \u2550\u2550\u2550\n';
    det+='Torque on pulley: T<sub>D</sub> = (T\u2081 \u2212 T\u2082) \u00D7 R<sub>D</sub>\n';
    det+='T\u2081/T\u2082 = '+tr+'\n';
    det+=T.toFixed(2)+' = ('+tr+' \u00D7 T\u2082 \u2212 T\u2082) \u00D7 ('+pDia+'/2)\n';
    det+='\u2234 T\u2082 = '+T2val.toFixed(2)+' N\n';
    det+='and T\u2081 = '+T1val.toFixed(2)+' N\n\n';
    det+='Thus vertical component = (T\u2081 + T\u2082) sin \u03B8 = ('+T1val.toFixed(2)+' + '+T2val.toFixed(2)+') \u00D7 sin('+thetaDeg+'\u00B0) = '+beltV.toFixed(2)+' N\n';
    det+='horizontal component = (T\u2081 + T\u2082) cos \u03B8 = ('+T1val.toFixed(2)+' + '+T2val.toFixed(2)+') \u00D7 cos('+thetaDeg+'\u00B0) = '+beltH.toFixed(2)+' N\n\n';

    det+='\u2550\u2550\u2550 Analysis of Gear C \u2550\u2550\u2550\n';
    det+='F<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = 2 \u00D7 '+T.toFixed(2)+' / '+gd+' = '+FtC.toFixed(2)+' N\n';
    det+='F<sub>rC</sub> = F<sub>tC</sub> tan \u03B1 = '+FtC.toFixed(2)+' \u00D7 tan '+alphaDeg+'\u00B0 = '+FrC.toFixed(2)+' N\n\n';

    // Vertical loading
    det+='\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\n';
    det+='Taking moments about bearing A, and equating to zero:\n';
    det+='R<sub>BV</sub>(AB) = [F<sub>tC</sub> + W<sub>C</sub>](AC) + [(T\u2081+T\u2082) sin \u03B8 + W<sub>D</sub>](AD)\n';
    det+='R<sub>BV</sub> \u00D7 '+AB+' = [('+FtC.toFixed(2)+' + '+Wc+') \u00D7 '+AC+'] + [('+beltV.toFixed(2)+' + '+Wd+') \u00D7 '+AD+']\n';
    det+='\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\n';
    det+='Also R<sub>AV</sub> + R<sub>BV</sub> = [F<sub>tC</sub> + W<sub>C</sub>] + [(T\u2081+T\u2082) sin \u03B8 + W<sub>D</sub>]\n';
    det+='R<sub>AV</sub> + '+R_BV.toFixed(2)+' = ('+Fv_C.toFixed(2)+') + ('+Fv_D.toFixed(2)+')\n';
    det+='\u2234 R<sub>AV</sub> = '+R_AV.toFixed(2)+' N\n\n';

    det+='\u2022 Moments:\n';
    det+='M<sub>AV</sub> = M<sub>DV</sub> = 0\n';
    det+='M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+R_AV.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    if(AD>AB){
        det+='M<sub>BV</sub> = \u2212[(T\u2081+T\u2082) sin \u03B8 + W<sub>D</sub>] \u00D7 '+(AD-AB)+' = \u2212('+Fv_D.toFixed(2)+') \u00D7 '+(AD-AB)+' = '+M_BV.toFixed(2)+' N\u00B7mm\n';
    }
    det+='\n';

    // Horizontal loading
    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\n';
    det+='Taking moments about bearing A, and equating to zero:\n';
    det+='R<sub>BH</sub>(AB) = F<sub>rC</sub>(AC) + [(T\u2081+T\u2082) cos \u03B8](AD)\n';
    det+='R<sub>BH</sub> \u00D7 '+AB+' = ('+FrC.toFixed(2)+' \u00D7 '+AC+') + ('+beltH.toFixed(2)+' \u00D7 '+AD+')\n';
    det+='\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\n';
    det+='Also, R<sub>AH</sub> + R<sub>BH</sub> = F<sub>rC</sub> + [(T\u2081+T\u2082) cos \u03B8]\n';
    det+='R<sub>AH</sub> + '+R_BH.toFixed(2)+' = '+FrC.toFixed(2)+' + '+beltH.toFixed(2)+'\n';
    det+='\u2234 R<sub>AH</sub> = '+R_AH.toFixed(2)+' N\n';
    if(R_AH_neg){
        det+='\nSince the reaction is negative, R<sub>AH</sub> acts in opposite direction.\n';
    }
    det+='\n';

    det+='\u2022 Moments:\n';
    det+='M<sub>AH</sub> = M<sub>DH</sub> = 0\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+R_AH.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    if(AD>AB){
        det+='M<sub>BH</sub> = \u2212[(T\u2081+T\u2082) cos \u03B8] \u00D7 '+(AD-AB)+' = \u2212'+beltH.toFixed(2)+' \u00D7 '+(AD-AB)+' = '+M_BH.toFixed(2)+' N\u00B7mm\n';
    }
    det+='\n';

    // Resultant
    det+='\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\n';
    det+='Resultant moment at C:\n';
    det+='M<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2)\n';
    det+='= '+M_C.toFixed(2)+' N\u00B7mm\n\n';

    if(AD>AB){
        det+='Resultant moment at B:\n';
        det+='M<sub>B</sub> = \u221A(M\u00B2<sub>BV</sub> + M\u00B2<sub>BH</sub>) = \u221A(('+M_BV.toFixed(2)+')\u00B2 + ('+M_BH.toFixed(2)+')\u00B2)\n';
        det+='= '+M_B.toFixed(2)+' N\u00B7mm\n\n';
    }

    det+='Thus maximum bending moment occurs at \''+maxLoc+'\', i.e.\n';
    det+='<strong>M<sub>'+maxLoc+'</sub> = M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';

    return{x:x,L:totalL,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText,problemType){
    document.getElementById('results').style.display='block';document.getElementById('diagrams').style.display='block';
    document.getElementById('torqueResult').textContent=T.toFixed(2)+' N\u00B7mm ('+(T/1000).toFixed(2)+' N\u00B7m)';
    document.getElementById('momentResult').textContent=maxBM.toFixed(2)+' N\u00B7mm';
    if(shaftType==='solid')document.getElementById('diameterResult').textContent='d = '+outerDia+' mm';
    else document.getElementById('diameterResult').textContent='d\u2092 = '+outerDia+' mm, d\u1D62 = '+innerDia.toFixed(2)+' mm';
    var useShearOnly=(problemType==='twopulleyvh'||problemType==='pulleygear'||problemType==='singlepulley');
    var h='<h3>Step-by-Step Solution</h3><p><strong>Torque:</strong> T = '+T.toFixed(2)+' N\u00B7mm</p>';
    h+='<p><strong>\u03C3<sub>max</sub> = '+sigmaMax.toFixed(2)+' MPa, \u03C4<sub>max</sub> = '+tauMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>C<sub>m</sub> = '+Cm.toFixed(2)+', C<sub>t</sub> = '+Ct.toFixed(2)+'</strong></p>';
    if(document.getElementById('hasKeyway').checked)h+='<p><strong>\u26A0 Keyway applied (25% reduction)</strong></p>';
    h+='<h3>Force Analysis</h3><div style="background:#f4f1ff;padding:18px 20px;border-radius:12px;line-height:2;font-size:0.93em;border:1px solid #e0d8f8;">'+result.details.replace(/\n/g,'<br>')+'</div>';
    h+='<h3>Maximum Bending Moment</h3><p><strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong></p>';
    h+='<h3 style="color:#c0392b;">Diameter Calculation</h3>';
    if(useShearOnly){
        h+='<p><strong>According to maximum shear stress theory: For a solid shaft</strong></p>';
        h+='<p>d = [16/(\u03C0\u03C4<sub>max</sub>) \u00D7 \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)]<sup>1/3</sup></p>';
        h+='<p>= [16/(\u03C0 \u00D7 '+tauMax.toFixed(2)+') \u00D7 \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(2)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(2)+')\u00B2)]<sup>1/3</sup></p>';
        h+='<p><strong>\u2234 d = '+dShear.toFixed(2)+' mm</strong></p>';
    }else{
        h+='<p><strong style="color:#2c3e50;">Eq.(i) Max Normal Stress Theory:</strong></p>';
        h+='<p>d = [16/(\u03C0\u03C3<sub>max\u22121</sub>) \u00D7 {C<sub>m</sub>M + \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)}]<sup>1/3</sup></p>';
        h+='<p>= [16/(\u03C0 \u00D7 '+sigmaMax.toFixed(2)+') \u00D7 {'+Cm+' \u00D7 '+maxBM.toFixed(2)+' + \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(2)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(2)+')\u00B2)}]<sup>1/3</sup></p>';
        h+='<p><strong style="color:#8e44ad;">\u2234 d = '+dNormal.toFixed(2)+' mm ...Eq.(i)</strong></p><br>';
        h+='<p><strong style="color:#2c3e50;">Eq.(ii) Max Shear Stress Theory:</strong></p>';
        h+='<p>d = [16/(\u03C0\u03C4<sub>max\u22121</sub>) \u00D7 \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)]<sup>1/3</sup></p>';
        h+='<p>= [16/(\u03C0 \u00D7 '+tauMax.toFixed(2)+') \u00D7 \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(2)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(2)+')\u00B2)]<sup>1/3</sup></p>';
        h+='<p><strong style="color:#8e44ad;">\u2234 d = '+dShear.toFixed(2)+' mm ...Eq.(ii)</strong></p><br>';
        h+='<p>Based on Eqs. (i) and (ii), select maximum diameter for design, i.e.</p>';
        h+='<p>d = '+Math.max(dNormal,dShear).toFixed(2)+' mm</p>';
    }
    h+='<p style="color:#c0392b;font-size:1.3em;font-weight:bold;margin-top:12px;padding:10px;background:#fef0f0;border-radius:8px;display:inline-block;">\u2234 Standard size of shaft, d = '+outerDia+' mm</p>';
    if(deflText)h+='<h3>Angular Deflection \u03B8</h3>'+deflText;
    document.getElementById('detailedResults').innerHTML=h;document.getElementById('results').scrollIntoView({behavior:'smooth'});
}

function generateDiagrams(result,T){
    var x=result.x;var n=x.length;var tq=[];for(var i=0;i<n;i++)tq.push(T);
    document.getElementById('diagrams').innerHTML='<h2>\u2550 Horizontal Loading Diagram (HLD)</h2><div id="sfd_h"></div><div id="bmd_h"></div><h2>\u2550 Vertical Loading Diagram (VLD)</h2><div id="sfd_v"></div><div id="bmd_v"></div><h2>\u2550 Resultant Bending Moment</h2><div id="bmd_res"></div><h2>\u2550 Torque Diagram</h2><div id="torque_d"></div>';
    var mH=0,iH=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_h[i])>mH){mH=Math.abs(result.BM_h[i]);iH=i;}}
    var mV=0,iV=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_v[i])>mV){mV=Math.abs(result.BM_v[i]);iV=i;}}
    var mR=0,iR=0;for(var i=0;i<n;i++){if(result.BM_res[i]>mR){mR=result.BM_res[i];iR=i;}}
    var ly={paper_bgcolor:'rgba(0,0,0,0)',plot_bgcolor:'#faf9ff',font:{family:'Segoe UI',size:12}};
    Plotly.newPlot('sfd_h',[{x:x,y:result.SF_h,type:'scatter',fill:'tozeroy',line:{color:'#3498db',width:2.5},fillcolor:'rgba(52,152,219,0.15)'}],Object.assign({},ly,{title:'Shear Force \u2014 Horizontal',xaxis:{title:'Distance (mm)'},yaxis:{title:'SF (N)'}}));
    Plotly.newPlot('bmd_h',[{x:x,y:result.BM_h,type:'scatter',fill:'tozeroy',line:{color:'#e74c3c',width:2.5},fillcolor:'rgba(231,76,60,0.15)'}],Object.assign({},ly,{title:'Bending Moment \u2014 Horizontal',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N\u00B7mm)'},annotations:[{x:x[iH],y:result.BM_h[iH],text:'M = '+result.BM_h[iH].toFixed(0),showarrow:true,arrowhead:2,ax:40,ay:-40,font:{color:'#e74c3c',size:11}}]}));
    Plotly.newPlot('sfd_v',[{x:x,y:result.SF_v,type:'scatter',fill:'tozeroy',line:{color:'#3498db',width:2.5},fillcolor:'rgba(52,152,219,0.15)'}],Object.assign({},ly,{title:'Shear Force \u2014 Vertical',xaxis:{title:'Distance (mm)'},yaxis:{title:'SF (N)'}}));
    Plotly.newPlot('bmd_v',[{x:x,y:result.BM_v,type:'scatter',fill:'tozeroy',line:{color:'#e74c3c',width:2.5},fillcolor:'rgba(231,76,60,0.15)'}],Object.assign({},ly,{title:'Bending Moment \u2014 Vertical',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N\u00B7mm)'},annotations:[{x:x[iV],y:result.BM_v[iV],text:'M = '+result.BM_v[iV].toFixed(0),showarrow:true,arrowhead:2,ax:40,ay:-40,font:{color:'#e74c3c',size:11}}]}));
    Plotly.newPlot('bmd_res',[{x:x,y:result.BM_res,type:'scatter',fill:'tozeroy',line:{color:'#9b59b6',width:2.5},fillcolor:'rgba(155,89,182,0.15)'}],Object.assign({},ly,{title:'Resultant BM = \u221A(M\u00B2h + M\u00B2v)',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N\u00B7mm)'},annotations:[{x:x[iR],y:mR,text:'M = '+mR.toFixed(0),showarrow:true,arrowhead:2,ax:40,ay:-40,font:{color:'#9b59b6',size:11}}]}));
    Plotly.newPlot('torque_d',[{x:x,y:tq,type:'scatter',fill:'tozeroy',line:{color:'#27ae60',width:2.5},fillcolor:'rgba(39,174,96,0.15)'}],Object.assign({},ly,{title:'Torque Diagram',xaxis:{title:'Distance (mm)'},yaxis:{title:'Torque (N\u00B7mm)'}}));
}
