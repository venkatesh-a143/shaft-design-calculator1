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
function toggleShaftType(){
    document.getElementById('hollowInputs').style.display=document.querySelector('input[name="shaftType"]:checked').value==='hollow'?'block':'none';
}
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
function toggleVHC(){
    var v=document.querySelector('input[name="vhC_tensionMethod"]:checked').value;
    document.getElementById('vhC_directInputs').style.display=v==='direct'?'block':'none';
    document.getElementById('vhC_frictionInputs').style.display=v==='friction'?'block':'none';
}
function toggleVHD(){
    var v=document.querySelector('input[name="vhD_tensionMethod"]:checked').value;
    document.getElementById('vhD_directInputs').style.display=v==='direct'?'block':'none';
    document.getElementById('vhD_frictionInputs').style.display=v==='friction'?'block':'none';
}
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
    if(useShearOnly){var term2_s=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));dShear=Math.pow((16/(Math.PI*tauMax))*term2_s,1/3);dNormal=dShear;}
    else{var term1=Cm*maxBM;var term2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));dNormal=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2),1/3);dShear=Math.pow((16/(Math.PI*tauMax))*term2,1/3);}
    var outerDia,innerDia;
    var stdSizes=[6,8,10,12,14,16,18,20,22,25,28,32,36,40,45,50,56,63,71,80,90,100,110,125,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,450,480,500,530,560,600];
    if(shaftType==='solid'){var dd=Math.max(dNormal,dShear);outerDia=null;for(var i=0;i<stdSizes.length;i++){if(stdSizes[i]>=dd){outerDia=stdSizes[i];break;}}if(outerDia===null)outerDia=Math.ceil(dd/10)*10;innerDia=0;}
    else{var k=evalFraction(document.getElementById('diameterRatio').value);if(useShearOnly){var term2_h=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));var dd2=Math.pow((16/(Math.PI*tauMax))*term2_h*(1/(1-Math.pow(k,4))),1/3);}else{var term1h=Cm*maxBM;var term2h=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));var doN=Math.pow((16/(Math.PI*sigmaMax))*(term1h+term2h)*(1/(1-Math.pow(k,4))),1/3);var doS=Math.pow((16/(Math.PI*tauMax))*term2h*(1/(1-Math.pow(k,4))),1/3);var dd2=Math.max(doN,doS);}outerDia=null;for(var i=0;i<stdSizes.length;i++){if(stdSizes[i]>=dd2){outerDia=stdSizes[i];break;}}if(outerDia===null)outerDia=Math.ceil(dd2/10)*10;innerDia=k*outerDia;}
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
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));var maxBM=0;for(var i=0;i<n;i++){if(BR[i]>maxBM)maxBM=BR[i];}
    var det='F<sub>t</sub> = 2T/D = 2 \u00D7 '+T.toFixed(2)+' / '+gd+' = '+Ft.toFixed(2)+' N\nF<sub>r</sub> = F<sub>t</sub> \u00D7 tan \u03B1 = '+Ft.toFixed(2)+' \u00D7 tan '+document.getElementById('normal_pressureAngle').value+'\u00B0 = '+Fr.toFixed(2)+' N\n\nHorizontal Reactions:\nR<sub>AH</sub> = '+RAh.toFixed(2)+' N\nR<sub>BH</sub> = '+RBh.toFixed(2)+' N\n\nVertical Reactions:\nR<sub>AV</sub> = '+RAv.toFixed(2)+' N\nR<sub>BV</sub> = '+RBv.toFixed(2)+' N';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcUDL(T){
    var L=parseFloat(document.getElementById('totalLength').value);var lU=parseFloat(document.getElementById('udlLength').value);var w=parseFloat(document.getElementById('udlIntensity').value);var b1=(L-lU)/2;var b2=b1+lU;var RA=(w*lU)/2;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);if(xi<b1||xi>b2){SFh.push(0);BMh.push(0);SFv.push(0);BMv.push(0);}else{var d=xi-b1;var sf=RA-w*d;var bm=RA*d-w*d*d/2;SFh.push(sf);BMh.push(bm);SFv.push(sf);BMv.push(bm);}}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));var maxBM=0;for(var i=0;i<n;i++){if(BR[i]>maxBM)maxBM=BR[i];}
    var det='w = '+w+' N/mm, l = '+lU+' mm\nW = w \u00D7 l = '+(w*lU).toFixed(2)+' N\nR<sub>A</sub> = R<sub>B</sub> = W/2 = '+RA.toFixed(2)+' N\nM<sub>max</sub> = wl\u00B2/8 = '+(w*lU*lU/8).toFixed(2)+' N\u00B7mm';
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcSinglePulley(T){var bc=document.querySelector('input[name="sp_bearingCount"]:checked').value;if(bc==='1')return calcSinglePulley1B(T);else return calcSinglePulley2B(T);}

function calcSinglePulley1B(T){
    var L=parseFloat(document.getElementById('sp_overhangDist').value);var pd=parseFloat(document.getElementById('sp_pulleyDia1').value);var Wp=parseFloat(document.getElementById('sp_pulleyWeight1').value);var R=pd/2;
    var T1,T2;var tM=document.querySelector('input[name="tensionInput"]:checked').value;var td='';
    if(tM==='direct'){T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);td='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';}
    else if(tM==='ratio'){var r=parseFloat(document.getElementById('sp_tensionRatio').value);T2=T/(R*(r-1));T1=r*T2;td='T\u2081/T\u2082 = '+r+'\nT = (T\u2081 \u2212 T\u2082) \u00D7 R\n'+T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\n(T\u2081 \u2212 T\u2082) = '+(T/R).toFixed(2)+' N\nAlso T\u2081 = '+r+' \u00D7 T\u2082\n('+r+' \u00D7 T\u2082 \u2212 T\u2082) = '+(T/R).toFixed(2)+'\nT\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+r+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';}
    else{var mu=parseFloat(document.getElementById('sp_mu').value);var tR=parseFloat(document.getElementById('sp_theta').value)*Math.PI/180;var ratio=Math.exp(mu*tR);T2=T/(R*(ratio-1));T1=ratio*T2;td='T = (T\u2081 \u2212 T\u2082) \u00D7 R\n'+T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\n(T\u2081 \u2212 T\u2082) = '+(T/R).toFixed(2)+' N\n\nAlso T\u2081/T\u2082 = e^(\u03BC\u03B8) = '+ratio.toFixed(4)+'\nT\u2081 = '+ratio.toFixed(4)+' \u00D7 T\u2082\n'+(ratio-1).toFixed(4)+' \u00D7 T\u2082 = '+(T/R).toFixed(2)+'\nT\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+ratio.toFixed(4)+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';}
    var tf=T1+T2+Wp;var maxBM=tf*L;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];for(var i=0;i<n;i++){var xi=(i/(n-1))*L;x.push(xi);SFh.push(0);BMh.push(0);SFv.push(tf);BMv.push(tf*xi);}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.abs(BMv[i]));
    var det='\u2550\u2550\u2550 Single Pulley \u2014 1 Bearing (Cantilever) \u2550\u2550\u2550\n\nTorque: T = '+T.toFixed(2)+' N\u00B7mm\n\n'+td+'\nMaximum bending moment:\nM = (T\u2081 + T\u2082 + W) \u00D7 L\n= ('+T1.toFixed(2)+' + '+T2.toFixed(2)+' + '+Wp+') \u00D7 '+L+'\n= '+tf.toFixed(2)+' \u00D7 '+L+'\n<strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcSinglePulley2B(T){
    var AB=parseFloat(document.getElementById('sp_bearingDist').value);var AC=parseFloat(document.getElementById('sp_pulleyPos').value);var pd=parseFloat(document.getElementById('sp_pulleyDia2').value);var Wp=parseFloat(document.getElementById('sp_pulleyWeight2').value);var R=pd/2;
    var T1,T2;var tM=document.querySelector('input[name="tensionInput"]:checked').value;var td='';
    if(tM==='direct'){T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);td='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';}
    else if(tM==='ratio'){var r=parseFloat(document.getElementById('sp_tensionRatio').value);T2=T/(R*(r-1));T1=r*T2;td='But torque on pulley, T = (T\u2081 \u2212 T\u2082) \u00D7 R\n'+T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\nT\u2081 \u2212 T\u2082 = '+(T/R).toFixed(2)+' N\n\nAlso T\u2081 = '+r+' \u00D7 T\u2082\n'+r+' \u00D7 T\u2082 \u2212 T\u2082 = '+(T/R).toFixed(2)+'\nT\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+r+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';}
    else{var mu=parseFloat(document.getElementById('sp_mu').value);var tR=parseFloat(document.getElementById('sp_theta').value)*Math.PI/180;var ratio=Math.exp(mu*tR);T2=T/(R*(ratio-1));T1=ratio*T2;td='T = (T\u2081 \u2212 T\u2082) \u00D7 R\n'+T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\nT\u2081 \u2212 T\u2082 = '+(T/R).toFixed(2)+' N\n\nAlso T\u2081/T\u2082 = e^(\u03BC\u03B8) = '+ratio.toFixed(4)+'\n'+(ratio-1).toFixed(4)+' \u00D7 T\u2082 = '+(T/R).toFixed(2)+'\nT\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+ratio.toFixed(4)+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';}
    var R_BV=(Wp*AC)/AB;var R_AV=Wp-R_BV;var Fh=T1+T2;var R_BH=(Fh*AC)/AB;var R_AH=Fh-R_BH;
    var M_CV=R_AV*AC;var M_CH=R_AH*AC;var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);var maxBM=M_C;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);if(xi<AC){SFv.push(R_AV);BMv.push(R_AV*xi);SFh.push(R_AH);BMh.push(R_AH*xi);}else{SFv.push(R_AV-Wp);BMv.push(R_AV*xi-Wp*(xi-AC));SFh.push(R_AH-Fh);BMh.push(R_AH*xi-Fh*(xi-AC));}}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 Single Pulley \u2014 2 Bearings \u2550\u2550\u2550\n\nTorque: T = '+T.toFixed(2)+' N\u00B7mm\n\n'+td+'\n\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nTaking moments about bearing A:\nR<sub>BV</sub> \u00D7 '+AB+' = W \u00D7 '+AC+' = '+Wp+' \u00D7 '+AC+'\n\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\nAlso R<sub>AV</sub> + R<sub>BV</sub> = W\n\u2234 R<sub>AV</sub> = '+R_AV.toFixed(2)+' N\n\nMoments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\nM<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+R_AV.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nTaking moments about bearing A:\nR<sub>BH</sub> \u00D7 '+AB+' = (T\u2081 + T\u2082) \u00D7 '+AC+'\n= ('+T1.toFixed(2)+' + '+T2.toFixed(2)+') \u00D7 '+AC+'\n\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\nAlso R<sub>AH</sub> + R<sub>BH</sub> = (T\u2081 + T\u2082)\n\u2234 R<sub>AH</sub> = '+R_AH.toFixed(2)+' N\n\nMoments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\nM<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+R_AH.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Resultant Moment at C \u2550\u2550\u2550\nM<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2)\n<strong>= '+M_C.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcTwoGear(T){
    var AB=parseFloat(document.getElementById('bearingDistance').value);
    var AC=parseFloat(document.getElementById('gearC_pos').value);
    var AD=parseFloat(document.getElementById('gearD_pos').value);
    var diaInput=document.querySelector('input[name="gearDiaInput"]:checked').value;
    var Dc,Dd,diaDetail='';
    if(diaInput==='direct'){
        Dc=parseFloat(document.getElementById('gearC_dia').value);
        Dd=parseFloat(document.getElementById('gearD_dia').value);
    }else{
        var mC=parseFloat(document.getElementById('gearC_module').value);
        var zC=parseFloat(document.getElementById('gearC_teeth').value);
        var mD=parseFloat(document.getElementById('gearD_module').value);
        var zD=parseFloat(document.getElementById('gearD_teeth').value);
        Dc=mC*zC;Dd=mD*zD;
        diaDetail='Gear C: D<sub>C</sub> = m \u00D7 z = '+mC+' \u00D7 '+zC+' = '+Dc+' mm\n';
        diaDetail+='Gear D: D<sub>D</sub> = m \u00D7 z = '+mD+' \u00D7 '+zD+' = '+Dd+' mm\n\n';
    }
    var Wc=parseFloat(document.getElementById('gearC_weight').value);
    var Wd=parseFloat(document.getElementById('gearD_weight').value);
    var alpha=parseFloat(document.getElementById('twoGear_pressureAngle').value)*Math.PI/180;
    var forceConfig=document.querySelector('input[name="twoGearForceConfig"]:checked').value;

    var FtC=(2*T)/Dc;var FrC=FtC*Math.tan(alpha);
    var FtD=(2*T)/Dd;var FrD=FtD*Math.tan(alpha);

    var V_RB,V_RA,H_RB,H_RA;
    var vForceC,vForceD,hForceC,hForceD;
    var vertDetailRB='',vertDetailRA='',horzDetailRB='',horzDetailRA='';
    var M_CV,M_DV,M_CH,M_DH;

    if(forceConfig==='crossed'){
        // CROSSED: Gear C: FtC horizontal, FrC vertical; Gear D: FtD vertical, FrD horizontal
        // Vertical plane: FtD down at D, FrC up at C (opposing) — no weights case
        // With weights: (FtD+Wd) at D, (FrC-Wc)... actually weights always act down
        // From textbook: R_BV*AB = FtD*AD - FrC*AC (without weights)
        // With weights: forces at C vertical = FrC downward? Let me re-read textbook...
        // Actually in the textbook Problem 21 (no weights, crossed):
        //   Vertical: R_BV*AB = FtD*AD - FrC*AC
        //   This means FtD acts one direction, FrC acts opposite
        //   R_AV = FtD - FrC - R_BV (net vertical = FtD - FrC)

        // But in textbook with same-direction+weights:
        //   Both act downward => R_BV*AB = (FtC+Wc)*AC + (FtD+Wd)*AD

        // For crossed without weights:
        vForceC=FrC; // FrC at C (acts upward in Fig, opposing FtD)
        vForceD=FtD; // FtD at D (acts downward)
        hForceC=FtC; // FtC at C (horizontal)
        hForceD=FrD; // FrD at D (horizontal)

        // Vertical: Taking moments about A
        // R_BV * AB = FtD * AD - FrC * AC
        V_RB=(FtD*AD-FrC*AC)/AB;
        V_RA=FtD-FrC-V_RB;

        vertDetailRB+='R<sub>BV</sub>(AB) = F<sub>tD</sub>(AD) \u2212 F<sub>rC</sub>(AC)\n';
        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = ('+FtD.toFixed(2)+' \u00D7 '+AD+') \u2212 ('+FrC.toFixed(2)+' \u00D7 '+AC+')\n';
        vertDetailRB+='\u2234 R<sub>BV</sub> = '+V_RB.toFixed(2)+' N\n\n';
        vertDetailRA+='Also R<sub>AV</sub> + R<sub>BV</sub> = F<sub>tD</sub> \u2212 F<sub>rC</sub>\n';
        vertDetailRA+='R<sub>AV</sub> + '+V_RB.toFixed(2)+' = '+FtD.toFixed(2)+' \u2212 '+FrC.toFixed(2)+'\n';
        vertDetailRA+='\u2234 R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';

        // Horizontal: Taking moments about A
        // R_BH*AB + FrD*AD = FtC*AC
        var H_RB_raw=(FtC*AC-FrD*AD)/AB;
        var H_RB_neg=H_RB_raw<0;var H_RB_abs=Math.abs(H_RB_raw);
        H_RB=H_RB_raw;

        horzDetailRB+='R<sub>BH</sub>(AB) + F<sub>rD</sub>(AD) = F<sub>tC</sub>(AC)\n';
        horzDetailRB+='R<sub>BH</sub> \u00D7 '+AB+' + ('+FrD.toFixed(2)+' \u00D7 '+AD+') = ('+FtC.toFixed(2)+' \u00D7 '+AC+')\n';
        horzDetailRB+='\u2234 R<sub>BH</sub> = '+H_RB_raw.toFixed(2)+' N';
        if(H_RB_neg) horzDetailRB+='\nSince the reaction is negative, R<sub>BH</sub> acts in opposite direction\nThus R<sub>BH</sub> = '+H_RB_abs.toFixed(2)+' N';
        horzDetailRB+='\n\n';

        // From revised Fig: R_AH + FrD = FtC + R_BH (if R_BH reversed)
        // General: R_AH = FtC - FrD - R_BH (if signs handled)
        // But from textbook: R_AH + FrD = FtC + |R_BH| (when R_BH is reversed)
        if(H_RB_neg){
            H_RA=FtC+H_RB_abs-FrD;
            horzDetailRA+='Also from revised diagram:\nR<sub>AH</sub> + F<sub>rD</sub> = F<sub>tC</sub> + R<sub>BH</sub>\n';
            horzDetailRA+='R<sub>AH</sub> + '+FrD.toFixed(2)+' = ('+FtC.toFixed(2)+' + '+H_RB_abs.toFixed(2)+')\n';
        }else{
            H_RA=FtC-FrD-H_RB;
            horzDetailRA+='Also R<sub>AH</sub> + R<sub>BH</sub> + F<sub>rD</sub> = F<sub>tC</sub>\n';
            horzDetailRA+='R<sub>AH</sub> + '+H_RB.toFixed(2)+' + '+FrD.toFixed(2)+' = '+FtC.toFixed(2)+'\n';
        }
        horzDetailRA+='\u2234 R<sub>AH</sub> = '+H_RA.toFixed(2)+' N\n\n';

        // Moments
        M_CV=V_RA*AC;
        M_DV=V_RB*(AB-AD);
        M_CH=H_RA*AC;
        M_DH=H_RB*(AB-AD); // negative if H_RB is negative

        // SFD/BMD
        var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
        // For SFD vertical: forces at A=R_AV(up), C=-FrC(up means opposing FtD, actually in diagram FrC acts upward at C)
        // Wait — in crossed config vertical: FtD at D acts DOWN, FrC at C acts... the equation says R_BV*AB = FtD*AD - FrC*AC
        // This means FrC creates moment in OPPOSITE direction to FtD. So FrC acts UPWARD.
        // For SFD: R_AV at A, then no force until D where FtD acts down, and at C... FrC acts up
        // But that means at C, the "vertical load" is FrC upward. So it SUBTRACTS from shear.
        // Actually, the equilibrium is: R_AV(up) - FrC(down at C??) + ... no.
        // Let me think again. R_BV*AB = FtD*AD - FrC*AC
        // If both FtD and FrC acted downward: R_BV*AB = FtD*AD + FrC*AC (positive contributions)
        // But we have a MINUS, so they oppose. Since FtD is at larger distance, and result is positive R_BV, FrC must act UPWARD.
        // So vertical SFD: R_AV(up at A), FrC(UP at C — but this means reaction-like... )
        // Actually in the FBD: R_AV and R_BV are reactions at bearings (upward).
        // External forces: FtD at D (downward) and FrC at C (upward? or downward?)
        // From the load diagram Fig 5.12(b): FrC is shown as a DOWNWARD arrow at C going upward? No...
        // Looking at Fig 5.12(b) vertical load diagram: FtD is shown downward at D, F_rC appears upward at C
        // Actually the arrow for F_rC in vertical diagram points UP
        // So vertical FBD: R_AV(up) + FrC(up at C) = FtD(down at D) + R_BV(... wait)
        // Actually: R_AV + R_BV = FtD - FrC (net downward force is FtD-FrC if FrC is up)
        // No: if FrC is UP, then: R_AV(up) + FrC(up) + R_BV(up) = FtD(down)
        // But R_AV + R_BV = FtD - FrC => This means FrC acts DOWN too but in opposite direction to FtD's moment
        // Hmm, let me just follow the math:
        // Sum of forces: R_AV + R_BV = FtD - FrC (this is the net, so FrC opposes FtD in force, meaning FrC is UP)
        // Moments about A: R_BV*AB = FtD*AD - FrC*AC

        // For SFD from left to right:
        // At A: +R_AV (upward)
        // At C: FrC acts upward = +FrC in shear
        // Wait no, if FrC is an external upward force, then passing through C, shear drops by FrC
        // Convention: positive SF = upward on left. Starting from left:
        // 0 to C: SF = R_AV
        // C to D: SF = R_AV + FrC (FrC is upward, so it adds? No... if FrC acts upward at C, then going past C from left, the net upward force increases)
        // Actually: Standard convention: cut section, sum forces to the LEFT.
        // Before C: SF = R_AV (upward)
        // After C (FrC acts upward at C): SF = R_AV + FrC... but wait, FrC is an external load acting upward, which is same direction as reaction. That's weird for a gear force.

        // Let me just check: R_AV = FtD - FrC - R_BV = 435.63 (from textbook)
        // FrC = 1034.50, FtD = 9474.21
        // At C: SF should jump. Before C: SF = 435.63
        // After C: if FrC acts UP: SF = 435.63 + 1034.50 = 1470.13? No that doesn't make sense
        // Actually, in the VLD Fig 5.12(b), the arrows show:
        //   R_AV upward at A, FrC downward at C? And FtD downward at D, R_BV upward at B
        //   Wait the diagram shows F_rC pointing UP at C on the VLD!

        // Let me try: FrC points UPWARD at C (opposing the downward loads)
        // Then equilibrium: R_AV (up) + FrC (up at C) + R_BV (up at B) = FtD (down at D)
        // R_AV + FrC + R_BV = FtD → R_AV + R_BV = FtD - FrC ✓ (matches)
        // Moment about A: R_BV*AB + FrC*AC = FtD*AD... but textbook says R_BV*AB = FtD*AD - FrC*AC
        // These are the same: R_BV*AB = FtD*AD - FrC*AC ✓

        // OK so FrC is UPWARD at C. For SFD:
        // Before C: SF = R_AV
        // Just after C: FrC acts upward, so net upward force to left increases: SF = R_AV - (-FrC) = R_AV + FrC
        // Wait... I need to be careful. If FrC acts UPWARD at C, and we're looking at forces to the LEFT of a cut:
        // Before C, only R_AV acts (upward), SF = R_AV
        // After C, forces to left = R_AV + FrC (both upward), SF = R_AV + FrC
        // After D, forces to left = R_AV + FrC - FtD, SF = R_AV + FrC - FtD
        // Check at B: SF should = -R_BV: R_AV + FrC - FtD = -(FtD - FrC - R_AV) + 2*FrC... 
        // Actually: R_AV + FrC - FtD = (FtD-FrC-R_BV) + FrC - FtD = -R_BV ✓

        for(var i=0;i<n;i++){
            var xi=(i/(n-1))*AB;x.push(xi);
            var sv,mv,sh,mh;
            if(xi<AC){
                sv=V_RA;mv=V_RA*xi;
            }else if(xi<AD){
                sv=V_RA+FrC;mv=V_RA*xi+FrC*(xi-AC);
            }else{
                sv=V_RA+FrC-FtD;mv=V_RA*xi+FrC*(xi-AC)-FtD*(xi-AD);
            }
            // Horizontal: FtC at C (let's say rightward), FrD at D (rightward too)
            // From equilibrium: R_AH + R_BH = FtC - FrD (if FtC and FrD oppose)
            // Actually: R_BH*AB + FrD*AD = FtC*AC, and looking at Fig 5.12(c):
            // HLD shows: R_AH leftward at A, FtC rightward at C, FrD rightward at D, R_BH at B
            // Actually in Fig 5.12(d) (revised): R_AH and R_BH both leftward, FtC rightward at C, FrD rightward at D
            // Wait, Fig 5.12(c) shows FtC at C pointing LEFT and FrD at D pointing RIGHT or...
            // From the equations: R_BH*AB + FrD*AD = FtC*AC → both FrD and FtC create positive moments
            // That means FtC creates CW moment and FrD creates CCW moment (opposing)? No...
            // FtC at distance AC creates moment FtC*AC (let's say CW)
            // FrD at distance AD creates moment FrD*AD (in same direction? No, it's on the other side)
            // R_BH*AB = FtC*AC - FrD*AD: so FtC*AC > FrD*AD in this case
            // They act in the SAME direction (both rightward) but FtC at AC creates more CW moment
            
            // For SFD horizontal: R_AH at A, then -FtC at C (FtC opposes R_AH), then +FrD at D
            // Wait: R_AH + R_BH = FtC - FrD (net horizontal if FtC and FrD act in same direction but FrD opposes the net)
            // Hmm this is getting confusing. Let me just use the actual computed values:
            // R_AH acts to the left (or right), and I'll track sign
            
            if(xi<AC){
                sh=H_RA;mh=H_RA*xi;
            }else if(xi<AD){
                sh=H_RA-FtC;mh=H_RA*xi-FtC*(xi-AC);
            }else{
                sh=H_RA-FtC+FrD;mh=H_RA*xi-FtC*(xi-AC)+FrD*(xi-AD);
            }
            SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);
        }

    }else{
        // SAME: Both gears FtC and FtD act vertically, FrC and FrD act horizontally (with weights)
        vForceC=FtC+Wc;vForceD=FtD+Wd;
        hForceC=FrC;hForceD=FrD;

        V_RB=(vForceC*AC+vForceD*AD)/AB;V_RA=vForceC+vForceD-V_RB;

        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = (F<sub>tC</sub>+W<sub>C</sub>)(AC) + (F<sub>tD</sub>+W<sub>D</sub>)(AD)\n';
        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = ('+vForceC.toFixed(2)+' \u00D7 '+AC+') + ('+vForceD.toFixed(2)+' \u00D7 '+AD+')\n';
        vertDetailRB+='\u2234 R<sub>BV</sub> = '+V_RB.toFixed(2)+' N\n\n';
        vertDetailRA+='R<sub>AV</sub> + R<sub>BV</sub> = (F<sub>tC</sub> + W<sub>C</sub>) + (F<sub>tD</sub> + W<sub>D</sub>)\n';
        vertDetailRA+='R<sub>AV</sub> + '+V_RB.toFixed(2)+' = '+vForceC.toFixed(2)+' + '+vForceD.toFixed(2)+'\n';
        vertDetailRA+='\u2234 R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';

        var H_RB_raw=(FrC*AC-FrD*AD)/AB;
        var H_RB_neg=H_RB_raw<0;var H_RB_abs=Math.abs(H_RB_raw);
        H_RB=H_RB_raw;
        if(H_RB_neg){H_RA=FrC+H_RB_abs-FrD;}else{H_RA=FrC+H_RB_raw-FrD;}

        horzDetailRB+='R<sub>BH</sub> \u00D7 '+AB+' + F<sub>rD</sub>(AD) = F<sub>rC</sub>(AC)\n';
        horzDetailRB+='R<sub>BH</sub> \u00D7 '+AB+' + ('+FrD.toFixed(2)+' \u00D7 '+AD+') = ('+FrC.toFixed(2)+' \u00D7 '+AC+')\n';
        horzDetailRB+='\u2234 R<sub>BH</sub> = '+H_RB_raw.toFixed(2)+' N';
        if(H_RB_neg) horzDetailRB+='\nSince negative, R<sub>BH</sub> acts in opposite direction\nThus R<sub>BH</sub> = '+H_RB_abs.toFixed(2)+' N';
        horzDetailRB+='\n\n';
        horzDetailRA+='R<sub>AH</sub> + F<sub>rD</sub> = F<sub>rC</sub> + R<sub>BH</sub>\nR<sub>AH</sub> + '+FrD.toFixed(2)+' = ('+FrC.toFixed(2)+' + '+H_RB_abs.toFixed(2)+')\n';
        horzDetailRA+='\u2234 R<sub>AH</sub> = '+H_RA.toFixed(2)+' N\n\n';

        M_CV=V_RA*AC;M_DV=V_RB*(AB-AD);
        if(H_RA<0){M_CH=-Math.abs(H_RA)*AC;}else{M_CH=H_RA*AC;}
        if(H_RB_neg){M_DH=-H_RB_abs*(AB-AD);}else{M_DH=H_RB_raw*(AB-AD);}

        var hForceA=(H_RA<0)?-Math.abs(H_RA):H_RA;
        var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
        for(var i=0;i<n;i++){
            var xi=(i/(n-1))*AB;x.push(xi);
            if(xi<AC){SFv.push(V_RA);BMv.push(V_RA*xi);SFh.push(hForceA);BMh.push(hForceA*xi);}
            else if(xi<AD){SFv.push(V_RA-vForceC);BMv.push(V_RA*xi-vForceC*(xi-AC));SFh.push(hForceA-FrC);BMh.push(hForceA*xi-FrC*(xi-AC));}
            else{SFv.push(V_RA-vForceC+vForceD);BMv.push(V_RA*xi-vForceC*(xi-AC)+vForceD*(xi-AD));SFh.push(hForceA-FrC+FrD);BMh.push(hForceA*xi-FrC*(xi-AC)+FrD*(xi-AD));}
        }
    }

    // Common: moments and resultant
    if(forceConfig==='crossed'){
        M_CV=V_RA*AC;M_DV=V_RB*(AB-AD);
        M_CH=H_RA*AC;M_DH=H_RB*(AB-AD);
    }
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);
    var M_D=Math.sqrt(M_DV*M_DV+M_DH*M_DH);
    var maxBM=Math.max(M_C,M_D);

    var n2=x.length;
    var BR=[];for(var i=0;i<n2;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));

    // Build details
    var det='';
    if(diaDetail) det+=diaDetail;
    det+='Gear C: F<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = 2 \u00D7 '+T.toFixed(2)+' / '+Dc+' = '+FtC.toFixed(2)+' N\n';
    det+='F<sub>rC</sub> = F<sub>tC</sub> \u00D7 tan '+document.getElementById('twoGear_pressureAngle').value+'\u00B0 = '+FtC.toFixed(2)+' \u00D7 tan '+document.getElementById('twoGear_pressureAngle').value+'\u00B0 = '+FrC.toFixed(2)+' N\n\n';
    det+='Gear D: F<sub>tD</sub> = 2M<sub>t</sub>/D<sub>D</sub> = 2 \u00D7 '+T.toFixed(2)+' / '+Dd+' = '+FtD.toFixed(2)+' N\n';
    det+='F<sub>rD</sub> = F<sub>tD</sub> \u00D7 tan '+document.getElementById('twoGear_pressureAngle').value+'\u00B0 = '+FtD.toFixed(2)+' \u00D7 tan '+document.getElementById('twoGear_pressureAngle').value+'\u00B0 = '+FrD.toFixed(2)+' N\n\n';

    det+='\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nTaking moments about A:\n';
    det+=vertDetailRB;det+=vertDetailRA;
    det+='Moments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\n';
    det+='M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+V_RA.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+V_RB.toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DV.toFixed(2)+' N\u00B7mm\n\n';

    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nTaking moments about A:\n';
    det+=horzDetailRB;det+=horzDetailRA;
    det+='Moments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+H_RA.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DH</sub> = '+((forceConfig==='crossed')?'R<sub>BH</sub>':'\u2212R<sub>BH</sub>')+' \u00D7 '+(AB-AD)+' = '+((forceConfig==='crossed')?'':'\\u2212')+Math.abs(H_RB).toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';

    det+='\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\n';
    det+='M<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2) = '+M_C.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>D</sub> = \u221A(M\u00B2<sub>DV</sub> + M\u00B2<sub>DH</sub>) = \u221A('+M_DV.toFixed(2)+'\u00B2 + ('+M_DH.toFixed(2)+')\u00B2) = '+M_D.toFixed(2)+' N\u00B7mm\n\n';
    det+='Thus maximum bending moment occurs at \''+(M_C>M_D?'C':'D')+'\', i.e.\n';
    det+='<strong>M<sub>'+(M_C>M_D?'C':'D')+'</sub> = M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';

    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcPulleyGear(T){
    var AB=parseFloat(document.getElementById('pg_bearingDistance').value);var AC=parseFloat(document.getElementById('pulley_pos').value);var pd=parseFloat(document.getElementById('pulley_dia').value);var Wa=parseFloat(document.getElementById('pulley_weight').value);var tr=parseFloat(document.getElementById('tensionRatio').value);var AD=parseFloat(document.getElementById('gear_pos').value);var gd=parseFloat(document.getElementById('gear_dia').value);var Wb=parseFloat(document.getElementById('gear_weight').value);var alpha=parseFloat(document.getElementById('pg_pressureAngle').value)*Math.PI/180;
    var Rp=pd/2;var T2val=T/(Rp*(tr-1));var T1val=tr*T2val;var FtD=(2*T)/gd;var FrD=FtD*Math.tan(alpha);
    var Fv_D=FtD+Wb;var R_BV=(Fv_D*AD)/AB;var R_AV=Fv_D-R_BV;var Fh_C=T1val+T2val+Wa;var R_BH=(FrD*AD-Fh_C*AC)/AB;var R_AH=FrD-Fh_C-R_BH;var R_AH_neg=R_AH<0;
    var M_CV=R_AV*AC;var M_DV=R_BV*(AB-AD);var M_CH=R_AH*AC;var M_DH=R_BH*(AB-AD);
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);var M_D=Math.sqrt(M_DV*M_DV+M_DH*M_DH);var maxBM=Math.max(M_C,M_D);var maxLoc=(M_C>=M_D)?'C':'D';
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=(i/(n-1))*AB;x.push(xi);var sv=0,mv=0,sh=0,mh=0;sv+=R_AV;mv+=R_AV*xi;if(xi>=AD){sv-=Fv_D;mv-=Fv_D*(xi-AD);}sh+=R_AH;mh+=R_AH*xi;if(xi>=AC){sh-=Fh_C;mh-=Fh_C*(xi-AC);}if(xi>=AD){sh+=FrD;mh+=FrD*(xi-AD);}SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='T\u2081 = '+T1val.toFixed(2)+' N, T\u2082 = '+T2val.toFixed(2)+' N\nPulley Force (T\u2081+T\u2082';if(Wa>0)det+='+W<sub>A</sub>';det+=') = '+Fh_C.toFixed(2)+' N\n\nF<sub>tD</sub> = 2M<sub>t</sub>/D<sub>B</sub> = '+FtD.toFixed(2)+' N\nF<sub>rD</sub> = F<sub>tD</sub> \u00D7 tan \u03B1 = '+FrD.toFixed(2)+' N\n';if(Wb>0)det+='Gear Force (F<sub>tD</sub>+W<sub>B</sub>) = '+Fv_D.toFixed(2)+' N\n';
    det+='\n\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nTaking moments about bearing A:\nR<sub>BV</sub>(AB) = F<sub>tD</sub>(AD)\nR<sub>BV</sub> \u00D7 '+AB+' = '+Fv_D.toFixed(2)+' \u00D7 '+AD+'\n\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\nAlso R<sub>AV</sub> + R<sub>BV</sub> = F<sub>tD</sub>\nR<sub>AV</sub> + '+R_BV.toFixed(2)+' = '+Fv_D.toFixed(2)+'\n\u2234 R<sub>AV</sub> = '+R_AV.toFixed(2)+' N\n\nMoments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\nM<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\nM<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+M_DV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nTaking moments about bearing A:\nR<sub>BH</sub>(AB) + [T\u2081+T\u2082](AC) = F<sub>rD</sub>(AD)\nR<sub>BH</sub> \u00D7 '+AB+' + [('+T1val.toFixed(2)+' + '+T2val.toFixed(2)+') \u00D7 '+AC+'] = ('+FrD.toFixed(2)+' \u00D7 '+AD+')\n\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\nAlso R<sub>AH</sub> + R<sub>BH</sub> + [T\u2081+T\u2082] = F<sub>rD</sub>\nR<sub>AH</sub> + '+R_BH.toFixed(2)+' + ('+T1val.toFixed(2)+' + '+T2val.toFixed(2)+') = '+FrD.toFixed(2)+'\n\u2234 R<sub>AH</sub> = '+R_AH.toFixed(2)+' N';if(R_AH_neg)det+='\nSince the reaction is negative, R<sub>AH</sub> acts in opposite direction';
    det+='\n\nMoments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\nM<sub>CH</sub> = '+R_AH.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\nM<sub>DH</sub> = '+R_BH.toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\nM<sub>C</sub> = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2) = '+M_C.toFixed(2)+' N\u00B7mm\n\nM<sub>D</sub> = \u221A('+M_DV.toFixed(2)+'\u00B2 + ('+M_DH.toFixed(2)+')\u00B2) = '+M_D.toFixed(2)+' N\u00B7mm\n\nThus maximum bending moment occurs at \''+maxLoc+'\', i.e.\n<strong>M<sub>'+maxLoc+'</sub> = M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det,torqueOverride:T};
}

function calcTwoPulleyVH(T_auto){
    var AB=parseFloat(document.getElementById('vh_bearingDist').value);var posC=parseFloat(document.getElementById('vh_pulleyC_pos').value);var diaC=parseFloat(document.getElementById('vh_pulleyC_dia').value);var posD=parseFloat(document.getElementById('vh_pulleyD_pos').value);var diaD=parseFloat(document.getElementById('vh_pulleyD_dia').value);var RC=diaC/2;var RD=diaD/2;
    var T3,T4,dMS='';var dM=document.querySelector('input[name="vhD_tensionMethod"]:checked').value;
    if(dM==='direct'){T3=parseFloat(document.getElementById('vh_T3D').value);T4=parseFloat(document.getElementById('vh_T4D').value);}else{var T3f=parseFloat(document.getElementById('vh_T3D_f').value);var muD=parseFloat(document.getElementById('vh_muD').value);var tDr=parseFloat(document.getElementById('vh_thetaD').value)*Math.PI/180;var rD=Math.exp(muD*tDr);T3=T3f;T4=T3/rD;dMS='T\u2083/T\u2084 = e^(\u03BC\u03B8) = '+rD.toFixed(4)+'\n'+T3.toFixed(2)+'/T\u2084 = '+rD.toFixed(4)+'\nT\u2084 = '+T4.toFixed(2)+' N\n';}
    var T_R=(T3-T4)*RD;var T1,T2,cMS='';var cM=document.querySelector('input[name="vhC_tensionMethod"]:checked').value;
    if(cM==='direct'){T1=parseFloat(document.getElementById('vh_T1C').value);T2=parseFloat(document.getElementById('vh_T2C').value);}else{var muC=parseFloat(document.getElementById('vh_muC').value);var tCr=parseFloat(document.getElementById('vh_thetaC').value)*Math.PI/180;var rC=Math.exp(muC*tCr);T2=T_R/(RC*(rC-1));T1=rC*T2;cMS='T\u2081/T\u2082 = e^(\u03BC\u03B8) = '+rC.toFixed(4)+'\nTorque = T<sub>R</sub> = '+T_R.toFixed(2)+' N\u00B7mm\n(T\u2081 \u2212 T\u2082) \u00D7 '+RC+' = '+T_R.toFixed(2)+'\nT\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+T1.toFixed(2)+' N\n';}
    var tM2=document.querySelector('input[name="vhTorqueMethod"]:checked').value;var T_used=(tM2==='pulley')?T_R:T_auto;
    var Fv=T1+T2;var RBV=(Fv*posC)/AB;var RAV=Fv-RBV;var Fh=T3+T4;var RBH=(Fh*posD)/AB;var RAH=Fh-RBH;var RAHn=RAH<0;
    var MCV=RAV*posC;var MCH=RAH*posC;var MBV=0;var MBH=(posD>AB)?-Fh*(posD-AB):0;
    var MC=Math.sqrt(MCV*MCV+MCH*MCH);var MB=Math.sqrt(MBV*MBV+MBH*MBH);var maxBM=Math.max(MC,MB);var mL=(MC>=MB)?'C':'B';
    var sX=Math.min(0,posC)-10;var eX=Math.max(AB,posD)+10;var tL=eX-sX;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=sX+(i/(n-1))*tL;x.push(xi);var sv=0,mv=0,sh=0,mh=0;if(xi>=0){sv+=RAV;mv+=RAV*xi;}if(xi>=posC){sv-=Fv;mv-=Fv*(xi-posC);}if(xi>=AB){sv+=RBV;mv+=RBV*(xi-AB);}if(xi>=0){sh+=RAH;mh+=RAH*xi;}if(xi>=AB){sh+=RBH;mh+=RBH*(xi-AB);}if(xi>=posD){sh-=Fh;mh-=Fh*(xi-posD);}SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 Right Pulley D \u2550\u2550\u2550\n';if(dM==='friction')det+=dMS;else det+='T\u2083 = '+T3.toFixed(2)+' N, T\u2084 = '+T4.toFixed(2)+' N\n';
    det+='\nT<sub>R</sub> = (T\u2083\u2212T\u2084)\u00D7R<sub>D</sub> = '+T_R.toFixed(2)+' N\u00B7mm\n\n\u2550\u2550\u2550 Left Pulley C \u2550\u2550\u2550\n';if(cM==='friction')det+=cMS;else det+='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';
    det+='\n\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nForce at C = '+Fv.toFixed(2)+' N\nR<sub>BV</sub> = '+RBV.toFixed(2)+' N, R<sub>AV</sub> = '+RAV.toFixed(2)+' N\n\n\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nForce at D = '+Fh.toFixed(2)+' N\nR<sub>BH</sub> = '+RBH.toFixed(2)+' N, R<sub>AH</sub> = '+RAH.toFixed(2)+' N';if(RAHn)det+=' (opposite)';
    det+='\n\nMoments:\nM<sub>CV</sub> = '+MCV.toFixed(2)+' N\u00B7mm\nM<sub>CH</sub> = '+MCH.toFixed(2)+' N\u00B7mm\n';if(posD>AB)det+='M<sub>BH</sub> = '+MBH.toFixed(2)+' N\u00B7mm\n';
    det+='\n\u2550\u2550\u2550 Resultant \u2550\u2550\u2550\nM<sub>C</sub> = '+MC.toFixed(2)+' N\u00B7mm\nM<sub>B</sub> = '+MB.toFixed(2)+' N\u00B7mm\n\n<strong>Max BM at '+mL+': M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:tL,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det,torqueOverride:T_used};
}

function calcGearAngledPulley(T){
    var AB=parseFloat(document.getElementById('gap_bearingDist').value);var AC=parseFloat(document.getElementById('gap_gearPos').value);var gd=parseFloat(document.getElementById('gap_gearDia').value);var alpha=parseFloat(document.getElementById('gap_pressureAngle').value)*Math.PI/180;var Wc=parseFloat(document.getElementById('gap_gearWeight').value);var AD=parseFloat(document.getElementById('gap_pulleyPos').value);var pD=parseFloat(document.getElementById('gap_pulleyDia').value);var tD=parseFloat(document.getElementById('gap_beltAngle').value);var tR=tD*Math.PI/180;var Wd=parseFloat(document.getElementById('gap_pulleyWeight').value);var tr=parseFloat(document.getElementById('gap_tensionRatio').value);
    var Rp=pD/2;var T2v=T/(Rp*(tr-1));var T1v=tr*T2v;var FtC=(2*T)/gd;var FrC=FtC*Math.tan(alpha);var bS=T1v+T2v;var bV=bS*Math.sin(tR);var bH=bS*Math.cos(tR);
    var FvC=FtC+Wc;var FvD=bV+Wd;var RBV=(FvC*AC+FvD*AD)/AB;var RAV=FvC+FvD-RBV;var RBH=(FrC*AC+bH*AD)/AB;var RAH=FrC+bH-RBH;var RAHn=RAH<0;
    var MCV=RAV*AC;var MCH=RAH*AC;var MBV,MBH,MDV,MDH;if(AD>AB){MBV=-FvD*(AD-AB);MBH=-bH*(AD-AB);MDV=0;MDH=0;}else{MBV=0;MBH=0;MDV=RBV*(AB-AD);MDH=RBH*(AB-AD);}
    var MC=Math.sqrt(MCV*MCV+MCH*MCH);var MB=Math.sqrt(MBV*MBV+MBH*MBH);var MD=Math.sqrt(MDV*MDV+MDH*MDH);var maxBM=Math.max(MC,MB,MD);var mL='C';if(MB>MC&&MB>=MD)mL='B';if(MD>MC&&MD>MB)mL='D';
    var sX=Math.min(0,AC)-10;var eX=Math.max(AB,AD)+10;var tL=eX-sX;var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=sX+(i/(n-1))*tL;x.push(xi);var sv=0,mv=0,sh=0,mh=0;if(xi>=0){sv+=RAV;mv+=RAV*xi;}if(xi>=AC){sv-=FvC;mv-=FvC*(xi-AC);}if(xi>=AB){sv+=RBV;mv+=RBV*(xi-AB);}if(xi>=AD){sv-=FvD;mv-=FvD*(xi-AD);}if(xi>=0){sh+=RAH;mh+=RAH*xi;}if(xi>=AC){sh-=FrC;mh-=FrC*(xi-AC);}if(xi>=AB){sh+=RBH;mh+=RBH*(xi-AB);}if(xi>=AD){sh-=bH;mh-=bH*(xi-AD);}SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 Pulley D \u2550\u2550\u2550\nT\u2081/T\u2082 = '+tr+', T\u2082 = '+T2v.toFixed(2)+' N, T\u2081 = '+T1v.toFixed(2)+' N\nVertical = '+bV.toFixed(2)+' N, Horizontal = '+bH.toFixed(2)+' N\n\n\u2550\u2550\u2550 Gear C \u2550\u2550\u2550\nF<sub>tC</sub> = '+FtC.toFixed(2)+' N, F<sub>rC</sub> = '+FrC.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Vertical \u2550\u2550\u2550\nR<sub>BV</sub> = '+RBV.toFixed(2)+' N, R<sub>AV</sub> = '+RAV.toFixed(2)+' N\nM<sub>CV</sub> = '+MCV.toFixed(2)+' N\u00B7mm\n';if(AD>AB)det+='M<sub>BV</sub> = '+MBV.toFixed(2)+' N\u00B7mm\n';else det+='M<sub>DV</sub> = '+MDV.toFixed(2)+' N\u00B7mm\n';
    det+='\n\u2550\u2550\u2550 Horizontal \u2550\u2550\u2550\nR<sub>BH</sub> = '+RBH.toFixed(2)+' N, R<sub>AH</sub> = '+RAH.toFixed(2)+' N';if(RAHn)det+=' (opposite)';det+='\nM<sub>CH</sub> = '+MCH.toFixed(2)+' N\u00B7mm\n';if(AD>AB)det+='M<sub>BH</sub> = '+MBH.toFixed(2)+' N\u00B7mm\n';else det+='M<sub>DH</sub> = '+MDH.toFixed(2)+' N\u00B7mm\n';
    det+='\n\u2550\u2550\u2550 Resultant \u2550\u2550\u2550\nM<sub>C</sub> = '+MC.toFixed(2)+' N\u00B7mm\n';if(AD>AB)det+='M<sub>B</sub> = '+MB.toFixed(2)+' N\u00B7mm\n';else det+='M<sub>D</sub> = '+MD.toFixed(2)+' N\u00B7mm\n';
    det+='\n<strong>Max BM at '+mL+': M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:tL,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText,problemType){
    document.getElementById('results').style.display='block';document.getElementById('diagrams').style.display='block';
    document.getElementById('torqueResult').textContent=T.toFixed(2)+' N\u00B7mm ('+(T/1000).toFixed(2)+' N\u00B7m)';
    document.getElementById('momentResult').textContent=maxBM.toFixed(2)+' N\u00B7mm';
    if(shaftType==='solid')document.getElementById('diameterResult').textContent='d = '+outerDia+' mm';
    else document.getElementById('diameterResult').textContent='d\u2092 = '+outerDia+' mm, d\u1D62 = '+innerDia.toFixed(2)+' mm';
    var useShearOnly=(problemType==='twopulleyvh'||problemType==='pulleygear'||problemType==='singlepulley');
    var h='<h3>Step-by-Step Solution</h3><p><strong>Torque:</strong> T = '+T.toFixed(2)+' N\u00B7mm</p>';
    h+='<p><strong>\u03C3<sub>max</sub> = '+sigmaMax.toFixed(2)+' MPa</strong></p><p><strong>\u03C4<sub>max</sub> = '+tauMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>C<sub>m</sub> = '+Cm.toFixed(2)+', C<sub>t</sub> = '+Ct.toFixed(2)+'</strong></p>';
    if(document.getElementById('hasKeyway').checked)h+='<p><strong>\u26A0 Keyway: stresses reduced by 25%</strong></p>';
    h+='<h3>Force Analysis</h3><div style="background:#f0f0f0;padding:15px;border-radius:5px;line-height:1.8;font-size:14px;">'+result.details.replace(/\n/g,'<br>')+'</div>';
    h+='<h3>Maximum Bending Moment</h3><p><strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong></p><h3>Diameter Calculation</h3>';
    if(useShearOnly){
        h+='<p><strong>According to maximum shear stress theory: For a solid shaft</strong></p>';
        h+='<p>d = [16/(\u03C0\u03C4<sub>max</sub>) \u00D7 \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)]<sup>1/3</sup></p>';
        h+='<p>= [16/(\u03C0 \u00D7 '+tauMax.toFixed(2)+') \u00D7 \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(2)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(2)+')\u00B2)]<sup>1/3</sup></p>';
        h+='<p><strong>\u2234 d = '+dShear.toFixed(2)+' mm</strong></p>';
    }else{
        h+='<p><strong>Eq.(i) Max Normal Stress Theory:</strong></p>';
        h+='<p>d = [16/(\u03C0\u03C3<sub>max</sub>) \u00D7 {C<sub>m</sub>M + \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)}]<sup>1/3</sup></p>';
        h+='<p>= [16/(\u03C0 \u00D7 '+sigmaMax.toFixed(2)+') \u00D7 {'+Cm+' \u00D7 '+maxBM.toFixed(2)+' + \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(2)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(2)+')\u00B2)}]<sup>1/3</sup></p>';
        h+='<p><strong>\u2234 d = '+dNormal.toFixed(2)+' mm ...Eq.(i)</strong></p>';
        h+='<br><p><strong>Eq.(ii) Max Shear Stress Theory:</strong></p>';
        h+='<p>d = [16/(\u03C0\u03C4<sub>max</sub>) \u00D7 \u221A((C<sub>m</sub>M)\u00B2 + (C<sub>t</sub>T)\u00B2)]<sup>1/3</sup></p>';
        h+='<p>= [16/(\u03C0 \u00D7 '+tauMax.toFixed(2)+') \u00D7 \u221A(('+Cm+' \u00D7 '+maxBM.toFixed(2)+')\u00B2 + ('+Ct+' \u00D7 '+T.toFixed(2)+')\u00B2)]<sup>1/3</sup></p>';
        h+='<p><strong>\u2234 d = '+dShear.toFixed(2)+' mm ...Eq.(ii)</strong></p>';
        h+='<br><p>Based on Eqs. (i) and (ii), select maximum diameter:</p><p>d = '+Math.max(dNormal,dShear).toFixed(2)+' mm</p>';
    }
    h+='<p style="color:red;font-size:1.3em;font-weight:bold;">\u2234 Standard size of shaft, d = '+outerDia+' mm</p>';
    if(deflText)h+='<h3>Angular Deflection \u03B8</h3>'+deflText;
    document.getElementById('detailedResults').innerHTML=h;document.getElementById('results').scrollIntoView({behavior:'smooth'});
}

function generateDiagrams(result,T){
    var x=result.x;var n=x.length;var tq=[];for(var i=0;i<n;i++)tq.push(T);
    document.getElementById('diagrams').innerHTML='<h2>\u2550 HLD</h2><div id="sfd_h"></div><div id="bmd_h"></div><h2>\u2550 VLD</h2><div id="sfd_v"></div><div id="bmd_v"></div><h2>\u2550 Resultant BM</h2><div id="bmd_res"></div><h2>\u2550 Torque</h2><div id="torque_d"></div>';
    var mH=0,iH=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_h[i])>mH){mH=Math.abs(result.BM_h[i]);iH=i;}}
    var mV=0,iV=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_v[i])>mV){mV=Math.abs(result.BM_v[i]);iV=i;}}
    var mR=0,iR=0;for(var i=0;i<n;i++){if(result.BM_res[i]>mR){mR=result.BM_res[i];iR=i;}}
    Plotly.newPlot('sfd_h',[{x:x,y:result.SF_h,type:'scatter',fill:'tozeroy',line:{color:'rgb(55,128,191)',width:3}}],{title:'SFD \u2014 Horizontal',xaxis:{title:'mm'},yaxis:{title:'N'}});
    Plotly.newPlot('bmd_h',[{x:x,y:result.BM_h,type:'scatter',fill:'tozeroy',line:{color:'rgb(219,64,82)',width:3}}],{title:'BMD \u2014 Horizontal',xaxis:{title:'mm'},yaxis:{title:'N\u00B7mm'},annotations:[{x:x[iH],y:result.BM_h[iH],text:'M = '+result.BM_h[iH].toFixed(0),showarrow:true,arrowhead:2,ax:40,ay:-40}]});
    Plotly.newPlot('sfd_v',[{x:x,y:result.SF_v,type:'scatter',fill:'tozeroy',line:{color:'rgb(55,128,191)',width:3}}],{title:'SFD \u2014 Vertical',xaxis:{title:'mm'},yaxis:{title:'N'}});
    Plotly.newPlot('bmd_v',[{x:x,y:result.BM_v,type:'scatter',fill:'tozeroy',line:{color:'rgb(219,64,82)',width:3}}],{title:'BMD \u2014 Vertical',xaxis:{title:'mm'},yaxis:{title:'N\u00B7mm'},annotations:[{x:x[iV],y:result.BM_v[iV],text:'M = '+result.BM_v[iV].toFixed(0),showarrow:true,arrowhead:2,ax:40,ay:-40}]});
    Plotly.newPlot('bmd_res',[{x:x,y:result.BM_res,type:'scatter',fill:'tozeroy',line:{color:'rgb(148,0,211)',width:3}}],{title:'Resultant BM',xaxis:{title:'mm'},yaxis:{title:'N\u00B7mm'},annotations:[{x:x[iR],y:mR,text:'M = '+mR.toFixed(0),showarrow:true,arrowhead:2,ax:40,ay:-40}]});
    Plotly.newPlot('torque_d',[{x:x,y:tq,type:'scatter',fill:'tozeroy',line:{color:'rgb(50,171,96)',width:3}}],{title:'Torque',xaxis:{title:'mm'},yaxis:{title:'N\u00B7mm'}});
}
