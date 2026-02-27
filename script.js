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

function evalFraction(s){
    s=s.trim();
    if(s.indexOf('/')!==-1){var p=s.split('/');return parseFloat(p[0])/parseFloat(p[1]);}
    return parseFloat(s);
}

function calculateShaft(){
    var problemType=document.querySelector('input[name="problemType"]:checked').value;
    var P=parseFloat(document.getElementById('power').value);
    var N=parseFloat(document.getElementById('speed').value);
    var T=(9.55e6*P)/N;

    var stressMethod=document.querySelector('input[name="stressMethod"]:checked').value;
    var sigmaMax,tauMax;

    if(stressMethod==='direct'){
        tauMax=parseFloat(document.getElementById('tauAllow').value);
        sigmaMax=2*tauMax;
    }else if(stressMethod==='material'){
        var sigUt=parseFloat(document.getElementById('sigmaUt').value);
        var sigYt=parseFloat(document.getElementById('sigmaYt').value);
        var nP=parseFloat(document.getElementById('nPrime').value);
        var sigU=sigUt/nP;var sigY=sigYt/nP;
        sigmaMax=Math.min(0.36*sigU,0.6*sigY);
        tauMax=Math.min(0.18*sigU,0.3*sigY);
    }else if(stressMethod==='working'){
        sigmaMax=parseFloat(document.getElementById('sigmaWorking').value);
        tauMax=parseFloat(document.getElementById('tauWorking').value);
    }else if(stressMethod==='yield'){
        var sYt=parseFloat(document.getElementById('sigmaYtYield').value);
        var fos=parseFloat(document.getElementById('fosYield').value);
        sigmaMax=sYt/fos;
        tauMax=sigmaMax/2;
    }

    if(document.getElementById('hasKeyway').checked){
        sigmaMax=0.75*sigmaMax;
        tauMax=0.75*tauMax;
    }

    var Cm=parseFloat(document.getElementById('cm').value);
    var Ct=parseFloat(document.getElementById('ct').value);

    var result;
    if(problemType==='normal') result=calcNormalShaft(T);
    else if(problemType==='udl') result=calcUDL(T);
    else if(problemType==='singlepulley') result=calcSinglePulley(T);
    else if(problemType==='twogear') result=calcTwoGear(T);
    else if(problemType==='pulleygear') result=calcPulleyGear(T);
    else if(problemType==='twopulleyvh') result=calcTwoPulleyVH(T);
    else if(problemType==='gearangledpulley') result=calcGearAngledPulley(T);

    if((problemType==='twopulleyvh'||problemType==='pulleygear') && result.torqueOverride!==undefined){
        T=result.torqueOverride;
    }

    var shaftType=document.querySelector('input[name="shaftType"]:checked').value;
    var maxBM=result.maxBM;
    var useShearOnly=(problemType==='twopulleyvh'||problemType==='pulleygear'||problemType==='singlepulley');

    var dNormal,dShear;
    if(useShearOnly){
        var term2_s=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));
        dShear=Math.pow((16/(Math.PI*tauMax))*term2_s,1/3);
        dNormal=dShear;
    }else{
        var term1=Cm*maxBM;
        var term2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));
        dNormal=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2),1/3);
        dShear=Math.pow((16/(Math.PI*tauMax))*term2,1/3);
    }

    var outerDia,innerDia;
    var stdSizes=[6,8,10,12,14,16,18,20,22,25,28,32,36,40,45,50,56,63,71,80,90,100,110,125,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,450,480,500,530,560,600];

    if(shaftType==='solid'){
        var dd=Math.max(dNormal,dShear);
        outerDia=null;
        for(var i=0;i<stdSizes.length;i++){if(stdSizes[i]>=dd){outerDia=stdSizes[i];break;}}
        if(outerDia===null) outerDia=Math.ceil(dd/10)*10;
        innerDia=0;
    }else{
        var k=evalFraction(document.getElementById('diameterRatio').value);
        if(useShearOnly){
            var term2_h=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));
            var dd2=Math.pow((16/(Math.PI*tauMax))*term2_h*(1/(1-Math.pow(k,4))),1/3);
        }else{
            var term1h=Cm*maxBM;
            var term2h=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));
            var doN=Math.pow((16/(Math.PI*sigmaMax))*(term1h+term2h)*(1/(1-Math.pow(k,4))),1/3);
            var doS=Math.pow((16/(Math.PI*tauMax))*term2h*(1/(1-Math.pow(k,4))),1/3);
            var dd2=Math.max(doN,doS);
        }
        outerDia=null;
        for(var i=0;i<stdSizes.length;i++){if(stdSizes[i]>=dd2){outerDia=stdSizes[i];break;}}
        if(outerDia===null) outerDia=Math.ceil(dd2/10)*10;
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

    displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText,problemType);
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

    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*AB;x.push(xi);
        if(xi<gp){SFh.push(RAh);BMh.push(RAh*xi);SFv.push(RAv);BMv.push(RAv*xi);}
        else{SFh.push(RAh-Ft);BMh.push(RAh*xi-Ft*(xi-gp));SFv.push(RAv-vf);BMv.push(RAv*xi-vf*(xi-gp));}
    }
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var maxBM=0;for(var i=0;i<n;i++){if(BR[i]>maxBM)maxBM=BR[i];}

    var det='F<sub>t</sub> = 2T/D = 2 \u00D7 '+T.toFixed(2)+' / '+gd+' = '+Ft.toFixed(2)+' N\n';
    det+='F<sub>r</sub> = F<sub>t</sub> \u00D7 tan \u03B1 = '+Ft.toFixed(2)+' \u00D7 tan '+document.getElementById('normal_pressureAngle').value+'\u00B0 = '+Fr.toFixed(2)+' N\n\n';
    det+='Horizontal Reactions:\nR<sub>AH</sub> = '+RAh.toFixed(2)+' N\nR<sub>BH</sub> = '+RBh.toFixed(2)+' N\n\n';
    det+='Vertical Reactions:\nR<sub>AV</sub> = '+RAv.toFixed(2)+' N\nR<sub>BV</sub> = '+RBv.toFixed(2)+' N';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcUDL(T){
    var L=parseFloat(document.getElementById('totalLength').value);
    var lU=parseFloat(document.getElementById('udlLength').value);
    var w=parseFloat(document.getElementById('udlIntensity').value);
    var b1=(L-lU)/2;var b2=b1+lU;var RA=(w*lU)/2;

    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);
        if(xi<b1||xi>b2){SFh.push(0);BMh.push(0);SFv.push(0);BMv.push(0);}
        else{var d=xi-b1;var sf=RA-w*d;var bm=RA*d-w*d*d/2;SFh.push(sf);BMh.push(bm);SFv.push(sf);BMv.push(bm);}
    }
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var maxBM=0;for(var i=0;i<n;i++){if(BR[i]>maxBM)maxBM=BR[i];}

    var det='w = '+w+' N/mm, l = '+lU+' mm\nW = w \u00D7 l = '+(w*lU).toFixed(2)+' N\n';
    det+='R<sub>A</sub> = R<sub>B</sub> = W/2 = '+RA.toFixed(2)+' N\nM<sub>max</sub> = wl\u00B2/8 = '+(w*lU*lU/8).toFixed(2)+' N\u00B7mm';
    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

// ========== SINGLE PULLEY ==========
function calcSinglePulley(T){
    var bearingCount=document.querySelector('input[name="sp_bearingCount"]:checked').value;
    if(bearingCount==='1') return calcSinglePulley1B(T);
    else return calcSinglePulley2B(T);
}

// 1 Bearing: Cantilever. M = (T1+T2+W)*L. Shear stress formula only.
function calcSinglePulley1B(T){
    var L=parseFloat(document.getElementById('sp_overhangDist').value);
    var pd=parseFloat(document.getElementById('sp_pulleyDia1').value);
    var Wp=parseFloat(document.getElementById('sp_pulleyWeight1').value);
    var R=pd/2;

    var T1,T2;
    var tMethod=document.querySelector('input[name="tensionInput"]:checked').value;
    var tensionDetail='';
    if(tMethod==='direct'){
        T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);
        tensionDetail='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';
    }else if(tMethod==='ratio'){
        var r=parseFloat(document.getElementById('sp_tensionRatio').value);
        // T = (T1-T2)*R, T1 = r*T2
        // T = (r*T2-T2)*R = T2*(r-1)*R
        T2=T/(R*(r-1));T1=r*T2;
        tensionDetail='T\u2081/T\u2082 = '+r+'\n';
        tensionDetail+='T = (T\u2081 \u2212 T\u2082) \u00D7 R\n';
        tensionDetail+=T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\n';
        tensionDetail+='(T\u2081 \u2212 T\u2082) = '+(T/R).toFixed(2)+' N\n';
        tensionDetail+='Also T\u2081 = '+r+' \u00D7 T\u2082\n';
        tensionDetail+='('+r+' \u00D7 T\u2082 \u2212 T\u2082) = '+(T/R).toFixed(2)+'\n';
        tensionDetail+='T\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+r+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';
    }else{
        var mu=parseFloat(document.getElementById('sp_mu').value);
        var thetaDeg=parseFloat(document.getElementById('sp_theta').value);
        var thetaRad=thetaDeg*Math.PI/180;
        var ratio=Math.exp(mu*thetaRad);
        // T1/T2 = ratio, (T1-T2)*R = T
        T2=T/(R*(ratio-1));T1=ratio*T2;
        tensionDetail+='T = (T\u2081 \u2212 T\u2082) \u00D7 R\n';
        tensionDetail+=T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\n';
        tensionDetail+='(T\u2081 \u2212 T\u2082) = '+(T/R).toFixed(2)+' N\n\n';
        tensionDetail+='Also T\u2081/T\u2082 = e^(\u03BC\u03B8) = e^('+mu+' \u00D7 '+thetaRad.toFixed(4)+') = '+ratio.toFixed(4)+'\n';
        tensionDetail+='T\u2081 = '+ratio.toFixed(4)+' \u00D7 T\u2082\n';
        tensionDetail+='('+ratio.toFixed(4)+' \u00D7 T\u2082 \u2212 T\u2082) = '+(T/R).toFixed(2)+'\n';
        tensionDetail+=''+(ratio-1).toFixed(4)+' \u00D7 T\u2082 = '+(T/R).toFixed(2)+'\n';
        tensionDetail+='T\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+ratio.toFixed(4)+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';
    }

    var totalForce=T1+T2+Wp;
    var maxBM=totalForce*L;

    // SFD/BMD for cantilever
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);
        SFh.push(0);BMh.push(0);
        SFv.push(totalForce);BMv.push(totalForce*xi);
    }
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.abs(BMv[i]));

    var det='\u2550\u2550\u2550 Single Pulley \u2014 1 Bearing (Cantilever) \u2550\u2550\u2550\n\n';
    det+='Torque: T = '+T.toFixed(2)+' N\u00B7mm\n\n';
    det+=tensionDetail+'\n';
    det+='Maximum bending moment:\n';
    det+='M = (T\u2081 + T\u2082 + W) \u00D7 L\n';
    det+='= ('+T1.toFixed(2)+' + '+T2.toFixed(2)+' + '+Wp+') \u00D7 '+L+'\n';
    det+='= '+totalForce.toFixed(2)+' \u00D7 '+L+'\n';
    det+='<strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';

    return{x:x,L:L,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

// 2 Bearings: Vertical = W only, Horizontal = (T1+T2) only. Shear stress formula only.
function calcSinglePulley2B(T){
    var AB=parseFloat(document.getElementById('sp_bearingDist').value);
    var AC=parseFloat(document.getElementById('sp_pulleyPos').value);
    var pd=parseFloat(document.getElementById('sp_pulleyDia2').value);
    var Wp=parseFloat(document.getElementById('sp_pulleyWeight2').value);
    var R=pd/2;

    var T1,T2;
    var tMethod=document.querySelector('input[name="tensionInput"]:checked').value;
    var tensionDetail='';
    if(tMethod==='direct'){
        T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);
        tensionDetail='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';
    }else if(tMethod==='ratio'){
        var r=parseFloat(document.getElementById('sp_tensionRatio').value);
        T2=T/(R*(r-1));T1=r*T2;
        tensionDetail+='But torque on pulley, T = (T\u2081 \u2212 T\u2082) \u00D7 R\n';
        tensionDetail+=T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\n';
        tensionDetail+='T\u2081 \u2212 T\u2082 = '+(T/R).toFixed(2)+' N\n\n';
        tensionDetail+='Also T\u2081 = '+r+' \u00D7 T\u2082\n';
        tensionDetail+=r+' \u00D7 T\u2082 \u2212 T\u2082 = '+(T/R).toFixed(2)+'\n';
        tensionDetail+='T\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+r+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';
    }else{
        var mu=parseFloat(document.getElementById('sp_mu').value);
        var thetaDeg=parseFloat(document.getElementById('sp_theta').value);
        var thetaRad=thetaDeg*Math.PI/180;
        var ratio=Math.exp(mu*thetaRad);
        T2=T/(R*(ratio-1));T1=ratio*T2;
        tensionDetail+='T = (T\u2081 \u2212 T\u2082) \u00D7 R\n';
        tensionDetail+=T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+R+'\n';
        tensionDetail+='T\u2081 \u2212 T\u2082 = '+(T/R).toFixed(2)+' N\n\n';
        tensionDetail+='Also T\u2081/T\u2082 = e^(\u03BC\u03B8) = e^('+mu+' \u00D7 '+thetaRad.toFixed(4)+') = '+ratio.toFixed(4)+'\n';
        tensionDetail+='T\u2081 = '+ratio.toFixed(4)+' \u00D7 T\u2082\n';
        tensionDetail+=(ratio-1).toFixed(4)+' \u00D7 T\u2082 = '+(T/R).toFixed(2)+'\n';
        tensionDetail+='T\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+ratio.toFixed(4)+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';
    }

    // VERTICAL: Only W acts downward at C
    var R_BV=(Wp*AC)/AB;
    var R_AV=Wp-R_BV;

    // HORIZONTAL: Only (T1+T2) acts at C
    var Fh=T1+T2;
    var R_BH=(Fh*AC)/AB;
    var R_AH=Fh-R_BH;

    // Moments at C
    var M_CV=R_AV*AC;
    var M_CH=R_AH*AC;

    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);
    var maxBM=M_C;

    // SFD/BMD
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*AB;x.push(xi);
        if(xi<AC){
            SFv.push(R_AV);BMv.push(R_AV*xi);
            SFh.push(R_AH);BMh.push(R_AH*xi);
        }else{
            SFv.push(R_AV-Wp);BMv.push(R_AV*xi-Wp*(xi-AC));
            SFh.push(R_AH-Fh);BMh.push(R_AH*xi-Fh*(xi-AC));
        }
    }
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));

    var det='\u2550\u2550\u2550 Single Pulley \u2014 2 Bearings \u2550\u2550\u2550\n\n';
    det+='Torque: T = '+T.toFixed(2)+' N\u00B7mm\n\n';
    det+=tensionDetail+'\n';

    det+='\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\n';
    det+='Taking moments about bearing A:\n';
    det+='R<sub>BV</sub> \u00D7 '+AB+' = W \u00D7 '+AC+'\n';
    det+='= '+Wp+' \u00D7 '+AC+'\n';
    det+='\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\n';
    det+='Also R<sub>AV</sub> + R<sub>BV</sub> = W\n';
    det+='R<sub>AV</sub> + '+R_BV.toFixed(2)+' = '+Wp+'\n';
    det+='\u2234 R<sub>AV</sub> = '+R_AV.toFixed(2)+' N\n\n';

    det+='Moments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\n';
    det+='M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+R_AV.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n\n';

    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\n';
    det+='Taking moments about bearing A:\n';
    det+='R<sub>BH</sub> \u00D7 '+AB+' = (T\u2081 + T\u2082) \u00D7 '+AC+'\n';
    det+='= ('+T1.toFixed(2)+' + '+T2.toFixed(2)+') \u00D7 '+AC+'\n';
    det+='\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\n';
    det+='Also R<sub>AH</sub> + R<sub>BH</sub> = (T\u2081 + T\u2082)\n';
    det+='R<sub>AH</sub> + '+R_BH.toFixed(2)+' = ('+T1.toFixed(2)+' + '+T2.toFixed(2)+')\n';
    det+='\u2234 R<sub>AH</sub> = '+R_AH.toFixed(2)+' N\n\n';

    det+='Moments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+R_AH.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n\n';

    det+='\u2550\u2550\u2550 Resultant Moment at C \u2550\u2550\u2550\n';
    det+='M<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2)\n';
    det+='<strong>= '+M_C.toFixed(2)+' N\u00B7mm</strong>';

    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

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

    var V_RB,V_RA,vertDetailRB='',vertDetailRA='',vForceC,vForceD;
    if(Wc>0||Wd>0){
        V_RB=((FtC+Wc)*AC+(FtD+Wd)*AD)/AB;V_RA=(FtC+Wc)+(FtD+Wd)-V_RB;
        vForceC=FtC+Wc;vForceD=FtD+Wd;
        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = (F<sub>tC</sub>+W<sub>C</sub>)(AC) + (F<sub>tD</sub>+W<sub>D</sub>)(AD)\n';
        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = ('+(FtC+Wc).toFixed(2)+' \u00D7 '+AC+') + ('+(FtD+Wd).toFixed(2)+' \u00D7 '+AD+')\n';
        vertDetailRB+='\u2234 R<sub>BV</sub> = '+V_RB.toFixed(2)+' N\n\n';
        vertDetailRA+='R<sub>AV</sub> + R<sub>BV</sub> = (F<sub>tC</sub> + W<sub>C</sub>) + (F<sub>tD</sub> + W<sub>D</sub>)\n';
        vertDetailRA+='R<sub>AV</sub> + '+V_RB.toFixed(2)+' = ('+(FtC+Wc).toFixed(2)+') + ('+(FtD+Wd).toFixed(2)+')\n';
        vertDetailRA+='\u2234 R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';
    }else{
        V_RB=(FtD*AD-FrC*AC)/AB;V_RA=FtD-FrC-V_RB;vForceC=FrC;vForceD=FtD;
        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = F<sub>tD</sub>(AD) \u2212 F<sub>rC</sub>(AC)\n';
        vertDetailRB+='R<sub>BV</sub> \u00D7 '+AB+' = ('+FtD.toFixed(2)+' \u00D7 '+AD+') \u2212 ('+FrC.toFixed(2)+' \u00D7 '+AC+')\n';
        vertDetailRB+='\u2234 R<sub>BV</sub> = '+V_RB.toFixed(2)+' N\n\n';
        vertDetailRA+='Also R<sub>AV</sub> + R<sub>BV</sub> = F<sub>tD</sub> \u2212 F<sub>rC</sub>\n';
        vertDetailRA+='\u2234 R<sub>AV</sub> = '+V_RA.toFixed(2)+' N\n\n';
    }

    var H_RB_raw=(FrC*AC-FrD*AD)/AB;var H_RB_neg=H_RB_raw<0;var H_RB_abs=Math.abs(H_RB_raw);
    var H_RA;if(H_RB_neg){H_RA=FrC+H_RB_abs-FrD;}else{H_RA=FrC+H_RB_raw-FrD;}
    var H_RA_neg=H_RA<0;var H_RA_abs=Math.abs(H_RA);

    var M_CV=V_RA*AC;var M_DV=V_RB*(AB-AD);
    var M_CH,M_DH;
    if(H_RA_neg){M_CH=-H_RA_abs*AC;}else{M_CH=H_RA*AC;}
    if(H_RB_neg){M_DH=-H_RB_abs*(AB-AD);}else{M_DH=H_RB_raw*(AB-AD);}
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);
    var M_D=Math.sqrt(M_DV*M_DV+M_DH*M_DH);
    var maxBM=Math.max(M_C,M_D);

    var hForceA=H_RA_neg?-H_RA_abs:H_RA;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*AB;x.push(xi);
        if(xi<AC){SFv.push(V_RA);BMv.push(V_RA*xi);SFh.push(hForceA);BMh.push(hForceA*xi);}
        else if(xi<AD){SFv.push(V_RA-vForceC);BMv.push(V_RA*xi-vForceC*(xi-AC));SFh.push(hForceA-FrC);BMh.push(hForceA*xi-FrC*(xi-AC));}
        else{SFv.push(V_RA-vForceC+vForceD);BMv.push(V_RA*xi-vForceC*(xi-AC)+vForceD*(xi-AD));SFh.push(hForceA-FrC+FrD);BMh.push(hForceA*xi-FrC*(xi-AC)+FrD*(xi-AD));}
    }
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));

    var det='Gear C: F<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = '+FtC.toFixed(2)+' N\nF<sub>rC</sub> = F<sub>tC</sub> \u00D7 tan \u03B1 = '+FrC.toFixed(2)+' N\n\n';
    det+='Gear D: F<sub>tD</sub> = 2M<sub>t</sub>/D<sub>D</sub> = '+FtD.toFixed(2)+' N\nF<sub>rD</sub> = F<sub>tD</sub> \u00D7 tan \u03B1 = '+FrD.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nTaking moments about A:\n';det+=vertDetailRB;det+=vertDetailRA;
    det+='Moments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\n';
    det+='M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+V_RA.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+V_RB.toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nTaking moments about A:\nR<sub>BH</sub> \u00D7 '+AB+' + F<sub>rD</sub>(AD) = F<sub>rC</sub>(AC)\n';
    det+='R<sub>BH</sub> \u00D7 '+AB+' + ('+FrD.toFixed(2)+' \u00D7 '+AD+') = ('+FrC.toFixed(2)+' \u00D7 '+AC+')\n';
    det+='\u2234 R<sub>BH</sub> = '+H_RB_raw.toFixed(2)+' N';
    if(H_RB_neg) det+='\nSince negative, R<sub>BH</sub> acts in opposite direction\nThus R<sub>BH</sub> = '+H_RB_abs.toFixed(2)+' N';
    det+='\n\nR<sub>AH</sub> + F<sub>rD</sub> = F<sub>rC</sub> + R<sub>BH</sub>\nR<sub>AH</sub> + '+FrD.toFixed(2)+' = ('+FrC.toFixed(2)+' + '+H_RB_abs.toFixed(2)+')\n';
    det+='\u2234 R<sub>AH</sub> = '+H_RA.toFixed(2)+' N';
    if(H_RA_neg) det+='\nSince negative, R<sub>AH</sub> acts in opposite direction\nThus R<sub>AH</sub> = '+H_RA_abs.toFixed(2)+' N';
    det+='\n\nMoments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+H_RA_abs.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DH</sub> = \u2212R<sub>BH</sub> \u00D7 '+(AB-AD)+' = \u2212'+H_RB_abs.toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\n';
    det+='M<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + '+M_CH.toFixed(2)+'\u00B2) = '+M_C.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>D</sub> = \u221A(M\u00B2<sub>DV</sub> + M\u00B2<sub>DH</sub>) = \u221A('+M_DV.toFixed(2)+'\u00B2 + ('+M_DH.toFixed(2)+')\u00B2) = '+M_D.toFixed(2)+' N\u00B7mm\n\n';
    det+='Thus maximum bending moment occurs at \''+(M_C>M_D?'C':'D')+'\', i.e.\n';
    det+='<strong>M<sub>'+(M_C>M_D?'C':'D')+'</sub> = M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function calcPulleyGear(T){
    var AB=parseFloat(document.getElementById('pg_bearingDistance').value);
    var AC=parseFloat(document.getElementById('pulley_pos').value);
    var pd=parseFloat(document.getElementById('pulley_dia').value);
    var Wa=parseFloat(document.getElementById('pulley_weight').value);
    var tr=parseFloat(document.getElementById('tensionRatio').value);
    var AD=parseFloat(document.getElementById('gear_pos').value);
    var gd=parseFloat(document.getElementById('gear_dia').value);
    var Wb=parseFloat(document.getElementById('gear_weight').value);
    var alpha=parseFloat(document.getElementById('pg_pressureAngle').value)*Math.PI/180;

    var Rp=pd/2;var T2val=T/(Rp*(tr-1));var T1val=tr*T2val;
    var FtD=(2*T)/gd;var FrD=FtD*Math.tan(alpha);
    var Fv_D=FtD+Wb;var R_BV=(Fv_D*AD)/AB;var R_AV=Fv_D-R_BV;
    var Fh_C=T1val+T2val+Wa;var R_BH=(FrD*AD-Fh_C*AC)/AB;var R_AH=FrD-Fh_C-R_BH;
    var R_AH_neg=R_AH<0;
    var M_CV=R_AV*AC;var M_DV=R_BV*(AB-AD);var M_CH=R_AH*AC;var M_DH=R_BH*(AB-AD);
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);var M_D=Math.sqrt(M_DV*M_DV+M_DH*M_DH);
    var maxBM=Math.max(M_C,M_D);var maxLoc=(M_C>=M_D)?'C':'D';

    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*AB;x.push(xi);var sv=0,mv=0,sh=0,mh=0;
        sv+=R_AV;mv+=R_AV*xi;if(xi>=AD){sv-=Fv_D;mv-=Fv_D*(xi-AD);}
        sh+=R_AH;mh+=R_AH*xi;if(xi>=AC){sh-=Fh_C;mh-=Fh_C*(xi-AC);}if(xi>=AD){sh+=FrD;mh+=FrD*(xi-AD);}
        SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);
    }
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));

    var det='T\u2081 = '+T1val.toFixed(2)+' N, T\u2082 = '+T2val.toFixed(2)+' N\n';
    det+='Pulley Force (T\u2081+T\u2082';if(Wa>0) det+='+W<sub>A</sub>';det+=') = '+Fh_C.toFixed(2)+' N\n\n';
    det+='F<sub>tD</sub> = 2M<sub>t</sub>/D<sub>B</sub> = '+FtD.toFixed(2)+' N\nF<sub>rD</sub> = F<sub>tD</sub> \u00D7 tan \u03B1 = '+FrD.toFixed(2)+' N\n';
    if(Wb>0) det+='Gear Force (F<sub>tD</sub>+W<sub>B</sub>) = '+Fv_D.toFixed(2)+' N\n';
    det+='\n\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nTaking moments about bearing A:\n';
    det+='R<sub>BV</sub>(AB) = F<sub>tD</sub>(AD)\nR<sub>BV</sub> \u00D7 '+AB+' = '+Fv_D.toFixed(2)+' \u00D7 '+AD+'\n';
    det+='\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\nAlso R<sub>AV</sub> + R<sub>BV</sub> = F<sub>tD</sub>\n';
    det+='R<sub>AV</sub> + '+R_BV.toFixed(2)+' = '+Fv_D.toFixed(2)+'\n\u2234 R<sub>AV</sub> = '+R_AV.toFixed(2)+' N\n\n';
    det+='Moments:\nM<sub>AV</sub> = M<sub>BV</sub> = 0\n';
    det+='M<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+R_AV.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+R_BV.toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DV.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nTaking moments about bearing A:\n';
    det+='R<sub>BH</sub>(AB) + [T\u2081+T\u2082](AC) = F<sub>rD</sub>(AD)\n';
    det+='R<sub>BH</sub> \u00D7 '+AB+' + [('+T1val.toFixed(2)+' + '+T2val.toFixed(2)+') \u00D7 '+AC+'] = ('+FrD.toFixed(2)+' \u00D7 '+AD+')\n';
    det+='\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\n';
    det+='Also R<sub>AH</sub> + R<sub>BH</sub> + [T\u2081+T\u2082] = F<sub>rD</sub>\n';
    det+='R<sub>AH</sub> + '+R_BH.toFixed(2)+' + ('+T1val.toFixed(2)+' + '+T2val.toFixed(2)+') = '+FrD.toFixed(2)+'\n';
    det+='\u2234 R<sub>AH</sub> = '+R_AH.toFixed(2)+' N';
    if(R_AH_neg) det+='\nSince the reaction is negative, R<sub>AH</sub> acts in opposite direction';
    det+='\n\nMoments:\nM<sub>AH</sub> = M<sub>BH</sub> = 0\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+R_AH.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>DH</sub> = R<sub>BH</sub> \u00D7 '+(AB-AD)+' = '+R_BH.toFixed(2)+' \u00D7 '+(AB-AD)+' = '+M_DH.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\n';
    det+='Resultant moment at C:\nM<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2)\n= '+M_C.toFixed(2)+' N\u00B7mm\n\n';
    det+='Resultant moment at D:\nM<sub>D</sub> = \u221A(M\u00B2<sub>DV</sub> + M\u00B2<sub>DH</sub>) = \u221A('+M_DV.toFixed(2)+'\u00B2 + ('+M_DH.toFixed(2)+')\u00B2)\n= '+M_D.toFixed(2)+' N\u00B7mm\n\n';
    det+='Thus maximum bending moment occurs at \''+maxLoc+'\', i.e.\n';
    det+='<strong>M<sub>'+maxLoc+'</sub> = M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:AB,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det,torqueOverride:T};
}

function calcTwoPulleyVH(T_auto){
    var AB=parseFloat(document.getElementById('vh_bearingDist').value);
    var posC=parseFloat(document.getElementById('vh_pulleyC_pos').value);
    var diaC=parseFloat(document.getElementById('vh_pulleyC_dia').value);
    var posD=parseFloat(document.getElementById('vh_pulleyD_pos').value);
    var diaD=parseFloat(document.getElementById('vh_pulleyD_dia').value);
    var RC=diaC/2;var RD=diaD/2;
    var T3,T4,dMethod_str='';
    var dMethod=document.querySelector('input[name="vhD_tensionMethod"]:checked').value;
    if(dMethod==='direct'){T3=parseFloat(document.getElementById('vh_T3D').value);T4=parseFloat(document.getElementById('vh_T4D').value);}
    else{var T3D_f=parseFloat(document.getElementById('vh_T3D_f').value);var muD=parseFloat(document.getElementById('vh_muD').value);var thetaD_rad=parseFloat(document.getElementById('vh_thetaD').value)*Math.PI/180;var ratioD_val=Math.exp(muD*thetaD_rad);T3=T3D_f;T4=T3/ratioD_val;dMethod_str='T\u2083/T\u2084 = e^(\u03BC\u03B8) = e^('+muD+' \u00D7 '+thetaD_rad.toFixed(4)+') = '+ratioD_val.toFixed(4)+'\n'+T3.toFixed(2)+'/T\u2084 = '+ratioD_val.toFixed(4)+'\nT\u2084 = '+T4.toFixed(2)+' N\n';}
    var T_R=(T3-T4)*RD;
    var T1,T2,cMethod_str='';
    var cMethod=document.querySelector('input[name="vhC_tensionMethod"]:checked').value;
    if(cMethod==='direct'){T1=parseFloat(document.getElementById('vh_T1C').value);T2=parseFloat(document.getElementById('vh_T2C').value);}
    else{var muC=parseFloat(document.getElementById('vh_muC').value);var thetaC_rad=parseFloat(document.getElementById('vh_thetaC').value)*Math.PI/180;var ratioC_val=Math.exp(muC*thetaC_rad);T2=T_R/(RC*(ratioC_val-1));T1=ratioC_val*T2;cMethod_str='T\u2081/T\u2082 = e^(\u03BC\u03B8) = e^('+muC+' \u00D7 '+thetaC_rad.toFixed(4)+') = '+ratioC_val.toFixed(4)+'\nTorque on left pulley = T<sub>R</sub> = '+T_R.toFixed(2)+' N\u00B7mm\n(T\u2081 \u2212 T\u2082) \u00D7 R<sub>C</sub> = T<sub>R</sub>\n('+ratioC_val.toFixed(4)+' \u00D7 T\u2082 \u2212 T\u2082) \u00D7 '+RC+' = '+T_R.toFixed(2)+'\nT\u2082 = '+T2.toFixed(2)+' N\nT\u2081 = '+ratioC_val.toFixed(4)+' \u00D7 '+T2.toFixed(2)+' = '+T1.toFixed(2)+' N\n';}
    var torqueMethod=document.querySelector('input[name="vhTorqueMethod"]:checked').value;
    var T_used=(torqueMethod==='pulley')?T_R:T_auto;
    var Fv_C=T1+T2;var R_BV=(Fv_C*posC)/AB;var R_AV=Fv_C-R_BV;
    var Fh_D=T3+T4;var R_BH=(Fh_D*posD)/AB;var R_AH=Fh_D-R_BH;var R_AH_neg=R_AH<0;
    var M_CV=R_AV*posC;var M_CH=R_AH*posC;var M_BV=0;var M_BH;
    if(posD>AB){M_BH=-Fh_D*(posD-AB);}else{M_BH=0;}
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);var M_B_res=Math.sqrt(M_BV*M_BV+M_BH*M_BH);
    var maxBM=Math.max(M_C,M_B_res);var maxLoc=(M_C>=M_B_res)?'C':'B';
    var startX=Math.min(0,posC)-10;var endX=Math.max(AB,posD)+10;var totalL=endX-startX;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=startX+(i/(n-1))*totalL;x.push(xi);var sv=0,mv=0,sh=0,mh=0;
        if(xi>=0){sv+=R_AV;mv+=R_AV*xi;}if(xi>=posC){sv-=Fv_C;mv-=Fv_C*(xi-posC);}if(xi>=AB){sv+=R_BV;mv+=R_BV*(xi-AB);}
        if(xi>=0){sh+=R_AH;mh+=R_AH*xi;}if(xi>=AB){sh+=R_BH;mh+=R_BH*(xi-AB);}if(xi>=posD){sh-=Fh_D;mh-=Fh_D*(xi-posD);}
        SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 Analysis of Right Pulley D \u2550\u2550\u2550\n';
    if(dMethod==='friction'){det+=dMethod_str;}else{det+='T\u2083 = '+T3.toFixed(2)+' N, T\u2084 = '+T4.toFixed(2)+' N\n';}
    det+='\nTorque on right pulley:\nT<sub>R</sub> = (T\u2083 \u2212 T\u2084) \u00D7 R<sub>D</sub> = ('+T3.toFixed(2)+' \u2212 '+T4.toFixed(2)+') \u00D7 '+RD+'\n= '+T_R.toFixed(2)+' N\u00B7mm\n\n';
    det+='\u2550\u2550\u2550 Analysis of Left Pulley C \u2550\u2550\u2550\n';
    if(cMethod==='friction'){det+=cMethod_str;}else{det+='T\u2081 = '+T1.toFixed(2)+' N, T\u2082 = '+T2.toFixed(2)+' N\n';}
    det+='\n\u2550\u2550\u2550 Vertical Loading (Pulley C \u2014 Belt Vertical) \u2550\u2550\u2550\nForce at C = T\u2081 + T\u2082 = '+T1.toFixed(2)+' + '+T2.toFixed(2)+' = '+Fv_C.toFixed(2)+' N\n\n';
    det+='Taking moments about A:\nR<sub>BV</sub> \u00D7 '+AB+' = '+Fv_C.toFixed(2)+' \u00D7 '+posC+'\n\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\u2234 R<sub>AV</sub> = '+Fv_C.toFixed(2)+' \u2212 '+R_BV.toFixed(2)+' = '+R_AV.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Horizontal Loading (Pulley D \u2014 Belt Horizontal) \u2550\u2550\u2550\nForce at D = T\u2083 + T\u2084 = '+T3.toFixed(2)+' + '+T4.toFixed(2)+' = '+Fh_D.toFixed(2)+' N\n\n';
    det+='Taking moments about A:\nR<sub>BH</sub> \u00D7 '+AB+' = '+Fh_D.toFixed(2)+' \u00D7 '+posD+'\n\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\u2234 R<sub>AH</sub> = '+Fh_D.toFixed(2)+' \u2212 '+R_BH.toFixed(2)+' = '+R_AH.toFixed(2)+' N';
    if(R_AH_neg) det+=' (acts opposite direction)';
    det+='\n\nMoments:\nM<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+posC+' = '+R_AV.toFixed(2)+' \u00D7 '+posC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    det+='M<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+posC+' = '+R_AH.toFixed(2)+' \u00D7 '+posC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    if(posD>AB){det+='M<sub>BH</sub> = \u2212F<sub>D</sub> \u00D7 ('+posD+' \u2212 '+AB+') = \u2212'+Fh_D.toFixed(2)+' \u00D7 '+(posD-AB)+' = '+M_BH.toFixed(2)+' N\u00B7mm\n';}
    det+='\n\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\nResultant moment at C:\nM<sub>C</sub> = \u221A(M\u00B2<sub>CV</sub> + M\u00B2<sub>CH</sub>) = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2)\n= '+M_C.toFixed(2)+' N\u00B7mm\n\n';
    det+='Resultant moment at B:\nM<sub>B</sub> = \u221A(M\u00B2<sub>BV</sub> + M\u00B2<sub>BH</sub>) = \u221A(0\u00B2 + ('+M_BH.toFixed(2)+')\u00B2)\n= '+M_B_res.toFixed(2)+' N\u00B7mm\n\n';
    det+='<strong>Maximum BM at '+maxLoc+': M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:totalL,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det,torqueOverride:T_used};
}

function calcGearAngledPulley(T){
    var AB=parseFloat(document.getElementById('gap_bearingDist').value);
    var AC=parseFloat(document.getElementById('gap_gearPos').value);
    var gd=parseFloat(document.getElementById('gap_gearDia').value);
    var alpha=parseFloat(document.getElementById('gap_pressureAngle').value)*Math.PI/180;
    var Wc=parseFloat(document.getElementById('gap_gearWeight').value);
    var AD=parseFloat(document.getElementById('gap_pulleyPos').value);
    var pDia=parseFloat(document.getElementById('gap_pulleyDia').value);
    var thetaDeg=parseFloat(document.getElementById('gap_beltAngle').value);
    var thetaRad=thetaDeg*Math.PI/180;
    var Wd=parseFloat(document.getElementById('gap_pulleyWeight').value);
    var tr=parseFloat(document.getElementById('gap_tensionRatio').value);
    var Rp=pDia/2;var T2val=T/(Rp*(tr-1));var T1val=tr*T2val;
    var FtC=(2*T)/gd;var FrC=FtC*Math.tan(alpha);
    var beltSum=T1val+T2val;var beltV=beltSum*Math.sin(thetaRad);var beltH=beltSum*Math.cos(thetaRad);
    var Fv_C=FtC+Wc;var Fv_D=beltV+Wd;
    var R_BV=(Fv_C*AC+Fv_D*AD)/AB;var R_AV=Fv_C+Fv_D-R_BV;
    var R_BH=(FrC*AC+beltH*AD)/AB;var R_AH=FrC+beltH-R_BH;
    var R_AH_neg=R_AH<0;
    var M_CV=R_AV*AC;var M_CH=R_AH*AC;
    var M_BV,M_BH,M_DV,M_DH;
    if(AD>AB){M_BV=-Fv_D*(AD-AB);M_BH=-beltH*(AD-AB);M_DV=0;M_DH=0;}
    else{M_BV=0;M_BH=0;M_DV=R_BV*(AB-AD);M_DH=R_BH*(AB-AD);}
    var M_C=Math.sqrt(M_CV*M_CV+M_CH*M_CH);var M_B_res=Math.sqrt(M_BV*M_BV+M_BH*M_BH);var M_D_res=Math.sqrt(M_DV*M_DV+M_DH*M_DH);
    var maxBM=Math.max(M_C,M_B_res,M_D_res);var maxLoc='C';
    if(M_B_res>M_C&&M_B_res>=M_D_res) maxLoc='B';if(M_D_res>M_C&&M_D_res>M_B_res) maxLoc='D';
    var startX=Math.min(0,AC)-10;var endX=Math.max(AB,AD)+10;var totalL=endX-startX;
    var n=500;var x=[],SFh=[],BMh=[],SFv=[],BMv=[];
    for(var i=0;i<n;i++){var xi=startX+(i/(n-1))*totalL;x.push(xi);var sv=0,mv=0,sh=0,mh=0;
        if(xi>=0){sv+=R_AV;mv+=R_AV*xi;}if(xi>=AC){sv-=Fv_C;mv-=Fv_C*(xi-AC);}if(xi>=AB){sv+=R_BV;mv+=R_BV*(xi-AB);}if(xi>=AD){sv-=Fv_D;mv-=Fv_D*(xi-AD);}
        if(xi>=0){sh+=R_AH;mh+=R_AH*xi;}if(xi>=AC){sh-=FrC;mh-=FrC*(xi-AC);}if(xi>=AB){sh+=R_BH;mh+=R_BH*(xi-AB);}if(xi>=AD){sh-=beltH;mh-=beltH*(xi-AD);}
        SFv.push(sv);BMv.push(mv);SFh.push(sh);BMh.push(mh);}
    var BR=[];for(var i=0;i<n;i++) BR.push(Math.sqrt(BMh[i]*BMh[i]+BMv[i]*BMv[i]));
    var det='\u2550\u2550\u2550 Analysis of Pulley D \u2550\u2550\u2550\nT\u2081/T\u2082 = '+tr+'\nT<sub>D</sub> = (T\u2081 \u2212 T\u2082) \u00D7 R<sub>D</sub> = T\n'+T.toFixed(2)+' = (T\u2081 \u2212 T\u2082) \u00D7 '+Rp+'\nT\u2082 = '+T2val.toFixed(2)+' N\nT\u2081 = '+T1val.toFixed(2)+' N\n\n';
    det+='Vertical component = (T\u2081 + T\u2082) \u00D7 sin '+thetaDeg+'\u00B0 = '+beltSum.toFixed(2)+' \u00D7 sin '+thetaDeg+'\u00B0 = '+beltV.toFixed(2)+' N\n';
    det+='Horizontal component = (T\u2081 + T\u2082) \u00D7 cos '+thetaDeg+'\u00B0 = '+beltSum.toFixed(2)+' \u00D7 cos '+thetaDeg+'\u00B0 = '+beltH.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Analysis of Gear C \u2550\u2550\u2550\nF<sub>tC</sub> = 2M<sub>t</sub>/D<sub>C</sub> = 2 \u00D7 '+T.toFixed(2)+' / '+gd+' = '+FtC.toFixed(2)+' N\n';
    det+='F<sub>rC</sub> = F<sub>tC</sub> \u00D7 tan '+document.getElementById('gap_pressureAngle').value+'\u00B0 = '+FrC.toFixed(2)+' N\n\n';
    det+='\u2550\u2550\u2550 Vertical Loading \u2550\u2550\u2550\nTaking moments about bearing A:\n';
    det+='R<sub>BV</sub>(AB) = [F<sub>tC</sub> + W<sub>C</sub>](AC) + [(T\u2081+T\u2082) sin \u03B8 + W<sub>D</sub>](AD)\n';
    det+='R<sub>BV</sub> \u00D7 '+AB+' = [('+(FtC).toFixed(2)+' + '+Wc+') \u00D7 '+AC+'] + [('+beltV.toFixed(2)+' + '+Wd+') \u00D7 '+AD+']\n';
    det+='\u2234 R<sub>BV</sub> = '+R_BV.toFixed(2)+' N\n\nAlso R<sub>AV</sub> + R<sub>BV</sub> = '+Fv_C.toFixed(2)+' + '+Fv_D.toFixed(2)+'\n\u2234 R<sub>AV</sub> = '+R_AV.toFixed(2)+' N\n\n';
    det+='Moments:\nM<sub>AV</sub> = M<sub>DV</sub> = 0\nM<sub>CV</sub> = R<sub>AV</sub> \u00D7 '+AC+' = '+R_AV.toFixed(2)+' \u00D7 '+AC+' = '+M_CV.toFixed(2)+' N\u00B7mm\n';
    if(AD>AB){det+='M<sub>BV</sub> = '+M_BV.toFixed(2)+' N\u00B7mm\n';}else{det+='M<sub>DV</sub> = R<sub>BV</sub> \u00D7 '+(AB-AD)+' = '+M_DV.toFixed(2)+' N\u00B7mm\n';}
    det+='\n\u2550\u2550\u2550 Horizontal Loading \u2550\u2550\u2550\nTaking moments about bearing A:\n';
    det+='R<sub>BH</sub>(AB) = F<sub>rC</sub>(AC) + [(T\u2081+T\u2082) cos \u03B8](AD)\nR<sub>BH</sub> \u00D7 '+AB+' = ('+FrC.toFixed(2)+' \u00D7 '+AC+') + ('+beltH.toFixed(2)+' \u00D7 '+AD+')\n\u2234 R<sub>BH</sub> = '+R_BH.toFixed(2)+' N\n\n';
    det+='Also R<sub>AH</sub> + R<sub>BH</sub> = F<sub>rC</sub> + [(T\u2081+T\u2082) cos \u03B8]\nR<sub>AH</sub> + '+R_BH.toFixed(2)+' = '+FrC.toFixed(2)+' + '+beltH.toFixed(2)+'\n\u2234 R<sub>AH</sub> = '+R_AH.toFixed(2)+' N';
    if(R_AH_neg) det+='\nSince the reaction is negative, R<sub>AH</sub> acts in opposite direction';
    det+='\n\nMoments:\nM<sub>AH</sub> = M<sub>DH</sub> = 0\nM<sub>CH</sub> = R<sub>AH</sub> \u00D7 '+AC+' = '+R_AH.toFixed(2)+' \u00D7 '+AC+' = '+M_CH.toFixed(2)+' N\u00B7mm\n';
    if(AD>AB){det+='M<sub>BH</sub> = '+M_BH.toFixed(2)+' N\u00B7mm\n';}else{det+='M<sub>DH</sub> = R<sub>BH</sub> \u00D7 '+(AB-AD)+' = '+M_DH.toFixed(2)+' N\u00B7mm\n';}
    det+='\n\u2550\u2550\u2550 Resultant Moments \u2550\u2550\u2550\nResultant moment at C:\nM<sub>C</sub> = \u221A('+M_CV.toFixed(2)+'\u00B2 + ('+M_CH.toFixed(2)+')\u00B2) = '+M_C.toFixed(2)+' N\u00B7mm\n\n';
    if(AD>AB){det+='Resultant moment at B:\nM<sub>B</sub> = \u221A(('+M_BV.toFixed(2)+')\u00B2 + ('+M_BH.toFixed(2)+')\u00B2) = '+M_B_res.toFixed(2)+' N\u00B7mm\n\n';}
    else{det+='Resultant moment at D:\nM<sub>D</sub> = \u221A(('+M_DV.toFixed(2)+')\u00B2 + ('+M_DH.toFixed(2)+')\u00B2) = '+M_D_res.toFixed(2)+' N\u00B7mm\n\n';}
    det+='Thus maximum bending moment occurs at \''+maxLoc+'\', i.e.\n<strong>M<sub>'+maxLoc+'</sub> = M = '+maxBM.toFixed(2)+' N\u00B7mm</strong>';
    return{x:x,L:totalL,SF_h:SFh,BM_h:BMh,SF_v:SFv,BM_v:BMv,BM_res:BR,maxBM:maxBM,details:det};
}

function displayResults(T,maxBM,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflText,problemType){
    document.getElementById('results').style.display='block';
    document.getElementById('diagrams').style.display='block';
    document.getElementById('torqueResult').textContent=T.toFixed(2)+' N\u00B7mm ('+(T/1000).toFixed(2)+' N\u00B7m)';
    document.getElementById('momentResult').textContent=maxBM.toFixed(2)+' N\u00B7mm';
    if(shaftType==='solid'){document.getElementById('diameterResult').textContent='d = '+outerDia+' mm';}
    else{document.getElementById('diameterResult').textContent='d\u2092 = '+outerDia+' mm, d\u1D62 = '+innerDia.toFixed(2)+' mm';}

    var useShearOnly=(problemType==='twopulleyvh'||problemType==='pulleygear'||problemType==='singlepulley');

    var h='<h3>Step-by-Step Solution</h3>';
    h+='<p><strong>Torque:</strong> T = '+T.toFixed(2)+' N\u00B7mm</p>';
    h+='<p><strong>\u03C3<sub>max</sub> = '+sigmaMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>\u03C4<sub>max</sub> = '+tauMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>C<sub>m</sub> = '+Cm.toFixed(2)+', C<sub>t</sub> = '+Ct.toFixed(2)+'</strong></p>';
    if(document.getElementById('hasKeyway').checked) h+='<p><strong>\u26A0 Keyway: stresses reduced by 25%</strong></p>';
    h+='<h3>Force Analysis</h3>';
    h+='<div style="background:#f0f0f0;padding:15px;border-radius:5px;line-height:1.8;font-size:14px;">'+result.details.replace(/\n/g,'<br>')+'</div>';
    h+='<h3>Maximum Bending Moment</h3>';
    h+='<p><strong>M = '+maxBM.toFixed(2)+' N\u00B7mm</strong></p>';
    h+='<h3>Diameter Calculation</h3>';

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
        h+='<br><p>Based on Eqs. (i) and (ii), select maximum diameter:</p>';
        h+='<p>d = '+Math.max(dNormal,dShear).toFixed(2)+' mm</p>';
    }

    h+='<p style="color:red;font-size:1.3em;font-weight:bold;">\u2234 Standard size of shaft, d = '+outerDia+' mm</p>';
    if(deflText){h+='<h3>Angular Deflection \u03B8</h3>'+deflText;}
    document.getElementById('detailedResults').innerHTML=h;
    document.getElementById('results').scrollIntoView({behavior:'smooth'});
}

function generateDiagrams(result,T){
    var x=result.x;var n=x.length;var tq=[];for(var i=0;i<n;i++)tq.push(T);
    document.getElementById('diagrams').innerHTML='<h2>\u2550 Horizontal Loading Diagram (HLD)</h2><div id="sfd_h"></div><div id="bmd_h"></div><h2>\u2550 Vertical Loading Diagram (VLD)</h2><div id="sfd_v"></div><div id="bmd_v"></div><h2>\u2550 Resultant Bending Moment</h2><div id="bmd_res"></div><h2>\u2550 Torque Diagram</h2><div id="torque_d"></div>';
    var mH=0,iH=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_h[i])>mH){mH=Math.abs(result.BM_h[i]);iH=i;}}
    var mV=0,iV=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_v[i])>mV){mV=Math.abs(result.BM_v[i]);iV=i;}}
    var mR=0,iR=0;for(var i=0;i<n;i++){if(result.BM_res[i]>mR){mR=result.BM_res[i];iR=i;}}
    Plotly.newPlot('sfd_h',[{x:x,y:result.SF_h,type:'scatter',fill:'tozeroy',line:{color:'rgb(55,128,191)',width:3}}],{title:'Shear Force Diagram \u2014 Horizontal Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'SF (N)'}});
    Plotly.newPlot('bmd_h',[{x:x,y:result.BM_h,type:'scatter',fill:'tozeroy',line:{color:'rgb(219,64,82)',width:3}}],{title:'Bending Moment Diagram \u2014 Horizontal Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N\u00B7mm)'},annotations:[{x:x[iH],y:result.BM_h[iH],text:'M<sub>max</sub> = '+result.BM_h[iH].toFixed(0)+' N\u00B7mm',showarrow:true,arrowhead:2,ax:40,ay:-40}]});
    Plotly.newPlot('sfd_v',[{x:x,y:result.SF_v,type:'scatter',fill:'tozeroy',line:{color:'rgb(55,128,191)',width:3}}],{title:'Shear Force Diagram \u2014 Vertical Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'SF (N)'}});
    Plotly.newPlot('bmd_v',[{x:x,y:result.BM_v,type:'scatter',fill:'tozeroy',line:{color:'rgb(219,64,82)',width:3}}],{title:'Bending Moment Diagram \u2014 Vertical Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N\u00B7mm)'},annotations:[{x:x[iV],y:result.BM_v[iV],text:'M<sub>max</sub> = '+result.BM_v[iV].toFixed(0)+' N\u00B7mm',showarrow:true,arrowhead:2,ax:40,ay:-40}]});
    Plotly.newPlot('bmd_res',[{x:x,y:result.BM_res,type:'scatter',fill:'tozeroy',line:{color:'rgb(148,0,211)',width:3}}],{title:'Resultant BM = \u221A(M\u00B2h + M\u00B2v)',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N\u00B7mm)'},annotations:[{x:x[iR],y:mR,text:'M<sub>max</sub> = '+mR.toFixed(0)+' N\u00B7mm',showarrow:true,arrowhead:2,ax:40,ay:-40}]});
    Plotly.newPlot('torque_d',[{x:x,y:tq,type:'scatter',fill:'tozeroy',line:{color:'rgb(50,171,96)',width:3}}],{title:'Torque Diagram',xaxis:{title:'Distance (mm)'},yaxis:{title:'Torque (N\u00B7mm)'}});
}
