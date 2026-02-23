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
function toggleShaftType(){
    document.getElementById('hollowInputs').style.display=document.querySelector('input[name="shaftType"]:checked').value==='hollow'?'block':'none';
}
function toggleTensionInput(){
    var v=document.querySelector('input[name="tensionInput"]:checked').value;
    document.getElementById('directTensionInputs').style.display=v==='direct'?'block':'none';
    document.getElementById('ratioTensionInputs').style.display=v==='ratio'?'block':'none';
}
function toggleModulusG(){document.getElementById('modulusGInput').style.display=document.getElementById('hasModulusG').checked?'block':'none';}
function toggleShaftLength(){document.getElementById('shaftLengthInput').style.display=document.getElementById('hasShaftLength').checked?'block':'none';}
function toggleAngleTwist(){document.getElementById('angleTwistInput').style.display=document.getElementById('hasAngleTwist').checked?'block':'none';}

document.getElementById('calcDeflection').addEventListener('change',function(){
    document.getElementById('deflectionInputs').style.display=this.checked?'block':'none';
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
        var sigU=sigUt/nP;
        var sigY=sigYt/nP;
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
    else if(problemType==='multipulley') result=calcMultiPulley(T);

    var shaftType=document.querySelector('input[name="shaftType"]:checked').value;
    var maxBM=result.maxBM;
    var term1=Cm*maxBM;
    var term2=Math.sqrt(Math.pow(Cm*maxBM,2)+Math.pow(Ct*T,2));

    var dNormal=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2),1/3);
    var dShear=Math.pow((16/(Math.PI*tauMax))*term2,1/3);

    var outerDia,innerDia,diameter;
    var stdSizes=[10,12,16,20,25,30,35,40,45,50,55,56,60,63,65,70,75,80,85,90,95,100,110,120];

    if(shaftType==='solid'){
        diameter=Math.max(dNormal,dShear);
        outerDia=stdSizes.find(function(s){return s>=diameter;})||Math.ceil(diameter/10)*10;
        innerDia=0;
    }else{
        var k=evalFraction(document.getElementById('diameterRatio').value);
        var doNorm=Math.pow((16/(Math.PI*sigmaMax))*(term1+term2)*(1/(1-Math.pow(k,4))),1/3);
        var doShear=Math.pow((16/(Math.PI*tauMax))*term2*(1/(1-Math.pow(k,4))),1/3);
        outerDia=Math.max(doNorm,doShear);
        outerDia=stdSizes.find(function(s){return s>=outerDia;})||Math.ceil(outerDia/10)*10;
        innerDia=k*outerDia;
        diameter=outerDia;
    }

    var deflectionText='';
    if(document.getElementById('calcDeflection').checked){
        var Gd=parseFloat(document.getElementById('deflG').value)*1000;
        var Ld=parseFloat(document.getElementById('deflL').value);
        var D4=Math.pow(outerDia,4);
        var theta=(584*T*Ld)/(Gd*D4);
        deflectionText='Angular Deflection theta = (584 x T x L) / (G x d^4) = (584 x '+T.toFixed(2)+' x '+Ld+') / ('+Gd+' x '+outerDia+'^4) = '+theta.toFixed(4)+' degrees';
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
    var Ft=(2*T)/gd;
    var Fr=Ft*Math.tan(alpha);
    var RAh=Ft*(AB-gp)/AB;
    var RBh=Ft-RAh;
    var vf=Fr+Wg;
    var RAv=vf*(AB-gp)/AB;
    var RBv=vf-RAv;
    var L=AB;var n=500;
    var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);
        if(xi<gp){SF_h.push(RAh);BM_h.push(RAh*xi);SF_v.push(RAv);BM_v.push(RAv*xi);}
        else{SF_h.push(RAh-Ft);BM_h.push(RAh*xi-Ft*(xi-gp));SF_v.push(RAv-vf);BM_v.push(RAv*xi-vf*(xi-gp));}
    }
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'normal',
        details:'Tangential Force Ft = 2T/D = '+Ft.toFixed(2)+' N\nRadial Force Fr = Ft x tan(alpha) = '+Fr.toFixed(2)+' N\n\nHorizontal Reactions:\nRA_h = '+RAh.toFixed(2)+' N\nRB_h = '+RBh.toFixed(2)+' N\n\nVertical Reactions:\nRA_v = '+RAv.toFixed(2)+' N\nRB_v = '+RBv.toFixed(2)+' N\n\nBM at gear (H) = '+BM_h[Math.round(gp/(L/(n-1)))].toFixed(2)+' N-mm\nBM at gear (V) = '+BM_v[Math.round(gp/(L/(n-1)))].toFixed(2)+' N-mm'};
}

function calcUDL(T){
    var L=parseFloat(document.getElementById('totalLength').value);
    var lU=parseFloat(document.getElementById('udlLength').value);
    var w=parseFloat(document.getElementById('udlIntensity').value);
    var b1=(L-lU)/2;var b2=b1+lU;var RA=(w*lU)/2;
    var n=500;var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);
        if(xi<b1||xi>b2){SF_h.push(0);BM_h.push(0);SF_v.push(0);BM_v.push(0);}
        else{var d=xi-b1;var sf=RA-w*d;var bm=RA*d-w*d*d/2;SF_h.push(sf);BM_h.push(bm);SF_v.push(sf);BM_v.push(bm);}
    }
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'udl',
        details:'UDL: w = '+w+' N/mm, l = '+lU+' mm\nTotal Load W = '+w*lU+' N\nRA = RB = '+(w*lU/2).toFixed(2)+' N\nMax BM at center = wl^2/8 = '+(w*lU*lU/8).toFixed(2)+' N-mm'};
}

function calcSinglePulley(T){
    var AB=parseFloat(document.getElementById('sp_bearingDist').value);
    var pp=parseFloat(document.getElementById('sp_pulleyPos').value);
    var pd=parseFloat(document.getElementById('sp_pulleyDia').value);
    var Wp=parseFloat(document.getElementById('sp_pulleyWeight').value);
    var T1,T2;
    if(document.querySelector('input[name="tensionInput"]:checked').value==='direct'){
        T1=parseFloat(document.getElementById('sp_T1').value);T2=parseFloat(document.getElementById('sp_T2').value);
    }else{
        var ratio=parseFloat(document.getElementById('sp_tensionRatio').value);
        T2=T/((pd/2)*(ratio-1));T1=ratio*T2;
    }
    var bf=T1+T2+Wp;
    var RAv=bf*(AB-pp)/AB;var RBv=bf-RAv;
    var L=AB;var n=500;
    var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);SF_h.push(0);BM_h.push(0);
        if(xi<pp){SF_v.push(RAv);BM_v.push(RAv*xi);}
        else{SF_v.push(RAv-bf);BM_v.push(RAv*xi-bf*(xi-pp));}
    }
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.abs(BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'singlepulley',
        details:'T1 = '+T1.toFixed(2)+' N, T2 = '+T2.toFixed(2)+' N\nTotal Belt Force (T1+T2) = '+(T1+T2).toFixed(2)+' N\nPulley Weight = '+Wp+' N\nTotal Vertical Force = '+bf.toFixed(2)+' N\n\nRA = '+RAv.toFixed(2)+' N\nRB = '+RBv.toFixed(2)+' N\nBending Moment M = (T1+T2)*L = '+(bf*pp*(AB-pp)/AB).toFixed(2)+' N-mm'};
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

    // Vertical: Radial forces act
    var V_RB=((FrD+Wd)*AD-(FrC+Wc)*AC)/AB;
    var V_RA=(FrD+Wd)-(FrC+Wc)+V_RB;

    // Horizontal: Tangential forces act
    var H_RB=(FtC*AC+FtD*AD)/AB;
    var H_RA=FtC+FtD-H_RB;

    var L=AB;var n=500;
    var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);
        if(xi<AC){SF_v.push(V_RA);BM_v.push(V_RA*xi);SF_h.push(H_RA);BM_h.push(H_RA*xi);}
        else if(xi<AD){SF_v.push(V_RA-(FrC+Wc));BM_v.push(V_RA*xi-(FrC+Wc)*(xi-AC));SF_h.push(H_RA-FtC);BM_h.push(H_RA*xi-FtC*(xi-AC));}
        else{SF_v.push(V_RA-(FrC+Wc)+(FrD+Wd));BM_v.push(V_RA*xi-(FrC+Wc)*(xi-AC)+(FrD+Wd)*(xi-AD));SF_h.push(H_RA-FtC+FtD);BM_h.push(H_RA*xi-FtC*(xi-AC)+FtD*(xi-AD));}
    }
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'twogear',
        details:'Gear C: FtC = '+FtC.toFixed(2)+' N, FrC = '+FrC.toFixed(2)+' N\nGear D: FtD = '+FtD.toFixed(2)+' N, FrD = '+FrD.toFixed(2)+' N\n\nVertical Reactions:\nRA_v = '+V_RA.toFixed(2)+' N, RB_v = '+V_RB.toFixed(2)+' N\n\nHorizontal Reactions:\nRA_h = '+H_RA.toFixed(2)+' N, RB_h = '+H_RB.toFixed(2)+' N\n\nMoments at C:\nM_cv = '+(V_RA*AC).toFixed(2)+' N-mm\nM_ch = '+(H_RA*AC).toFixed(2)+' N-mm\nResultant M_c = '+Math.sqrt(Math.pow(V_RA*AC,2)+Math.pow(H_RA*AC,2)).toFixed(2)+' N-mm\n\nMoments at D:\nM_dv = '+(V_RA*AD-(FrC+Wc)*(AD-AC)).toFixed(2)+' N-mm\nM_dh = '+(H_RA*AD-FtC*(AD-AC)).toFixed(2)+' N-mm\nResultant M_d = '+Math.sqrt(Math.pow(V_RA*AD-(FrC+Wc)*(AD-AC),2)+Math.pow(H_RA*AD-FtC*(AD-AC),2)).toFixed(2)+' N-mm'};
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
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);
        if(xi<bC){SF_v.push(-vP);BM_v.push(-vP*xi);SF_h.push(0);BM_h.push(0);}
        else if(xi<bD){SF_v.push(-vP+RCv);BM_v.push(-vP*xi+RCv*(xi-bC));SF_h.push(RCh);BM_h.push(RCh*(xi-bC));}
        else if(xi<=gp){SF_v.push(-vP+RCv+RDv);BM_v.push(-vP*xi+RCv*(xi-bC)+RDv*(xi-bD));SF_h.push(RCh+RDh);BM_h.push(RCh*(xi-bC)+RDh*(xi-bD));}
        else{SF_v.push(0);BM_v.push(0);SF_h.push(0);BM_h.push(0);}
    }
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.sqrt(BM_h[i]*BM_h[i]+BM_v[i]*BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'pulleygear',
        details:'Pulley: T1 = '+T1.toFixed(2)+' N, T2 = '+T2.toFixed(2)+' N\nGear: Ft = '+FtB.toFixed(2)+' N, Fr = '+FrB.toFixed(2)+' N\n\nVertical: RC = '+RCv.toFixed(2)+' N, RD = '+RDv.toFixed(2)+' N\nHorizontal: RC = '+RCh.toFixed(2)+' N, RD = '+RDh.toFixed(2)+' N'};
}

function calcMultiPulley(T){
    var AB=parseFloat(document.getElementById('mp_bearingDist').value);
    var pA=parseFloat(document.getElementById('mp_pulleyA_pos').value);
    var dA=parseFloat(document.getElementById('mp_pulleyA_dia').value);
    var WpA=parseFloat(document.getElementById('mp_pulleyA_weight').value);
    var T1A=parseFloat(document.getElementById('mp_T1A').value);var T2A=parseFloat(document.getElementById('mp_T2A').value);
    var pB=parseFloat(document.getElementById('mp_pulleyB_pos').value);
    var dB=parseFloat(document.getElementById('mp_pulleyB_dia').value);
    var WpB=parseFloat(document.getElementById('mp_pulleyB_weight').value);
    var T1B=parseFloat(document.getElementById('mp_T1B').value);var T2B=parseFloat(document.getElementById('mp_T2B').value);

    var fA=T1A+T2A+WpA;var fB=T1B+T2B+WpB;
    var R2v=(fA*pA+fB*pB)/AB;var R1v=fA+fB-R2v;
    var L=AB;var n=500;
    var x=[],SF_h=[],BM_h=[],SF_v=[],BM_v=[];
    for(var i=0;i<n;i++){
        var xi=(i/(n-1))*L;x.push(xi);SF_h.push(0);BM_h.push(0);
        if(xi<pA){SF_v.push(R1v);BM_v.push(R1v*xi);}
        else if(xi<pB){SF_v.push(R1v-fA);BM_v.push(R1v*xi-fA*(xi-pA));}
        else{SF_v.push(R1v-fA-fB);BM_v.push(R1v*xi-fA*(xi-pA)-fB*(xi-pB));}
    }
    var BM_res=[];for(var i=0;i<n;i++)BM_res.push(Math.abs(BM_v[i]));
    var maxBM=Math.max.apply(null,BM_res);
    return{x:x,L:L,SF_h:SF_h,BM_h:BM_h,SF_v:SF_v,BM_v:BM_v,BM_res:BM_res,maxBM:maxBM,type:'multipulley',
        details:'Pulley A: T1='+T1A+'N, T2='+T2A+'N, W='+WpA+'N, Total='+fA.toFixed(2)+'N\nPulley B: T1='+T1B+'N, T2='+T2B+'N, W='+WpB+'N, Total='+fB.toFixed(2)+'N\n\nR1 = '+R1v.toFixed(2)+' N\nR2 = '+R2v.toFixed(2)+' N'};
}

function displayResults(T,maxBM,diameter,outerDia,innerDia,shaftType,Cm,Ct,sigmaMax,tauMax,dNormal,dShear,result,deflectionText){
    document.getElementById('results').style.display='block';
    document.getElementById('diagrams').style.display='block';
    document.getElementById('torqueResult').textContent=T.toFixed(2)+' N-mm ('+(T/1000).toFixed(2)+' N-m)';
    document.getElementById('momentResult').textContent=maxBM.toFixed(2)+' N-mm';
    if(shaftType==='solid') document.getElementById('diameterResult').textContent='d = '+outerDia+' mm';
    else document.getElementById('diameterResult').textContent='do = '+outerDia+' mm, di = '+innerDia.toFixed(2)+' mm';

    var h='<h3>Step-by-Step Solution</h3>';
    h+='<p><strong>Torque:</strong> T = (9.55 x 10^6 x P) / n = (9.55 x 10^6 x '+document.getElementById('power').value+') / '+document.getElementById('speed').value+' = '+T.toFixed(2)+' N-mm</p>';
    h+='<p><strong>sigma_max = '+sigmaMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>tau_max = '+tauMax.toFixed(2)+' MPa</strong></p>';
    h+='<p><strong>Cm = '+Cm.toFixed(2)+', Ct = '+Ct.toFixed(2)+'</strong></p>';
    h+='<h3>Force Analysis</h3>';
    h+='<pre style="background:#f0f0f0;padding:15px;border-radius:5px;white-space:pre-wrap;font-size:14px;">'+result.details+'</pre>';
    h+='<h3>Maximum Bending Moment</h3>';
    h+='<p><strong>M = '+maxBM.toFixed(2)+' N-mm</strong></p>';
    h+='<h3>Diameter Calculation</h3>';
    h+='<p><strong>Eq.(i) Max Normal Stress Theory:</strong></p>';
    h+='<p>d = [16/(pi x '+sigmaMax.toFixed(2)+') x {'+Cm+' x '+maxBM.toFixed(0)+' + sqrt(('+Cm+' x '+maxBM.toFixed(0)+')^2 + ('+Ct+' x '+T.toFixed(0)+')^2)}]^(1/3)</p>';
    h+='<p><strong>d = '+dNormal.toFixed(2)+' mm ...Eq.(i)</strong></p>';
    h+='<p><strong>Eq.(ii) Max Shear Stress Theory:</strong></p>';
    h+='<p>d = [16/(pi x '+tauMax.toFixed(2)+') x sqrt(('+Cm+' x '+maxBM.toFixed(0)+')^2 + ('+Ct+' x '+T.toFixed(0)+')^2)]^(1/3)</p>';
    h+='<p><strong>d = '+dShear.toFixed(2)+' mm ...Eq.(ii)</strong></p>';
    h+='<p>Select maximum diameter: d = '+Math.max(dNormal,dShear).toFixed(2)+' mm</p>';
    h+='<p style="color:red;font-size:1.3em;font-weight:bold;">Standard size of shaft, d = '+outerDia+' mm</p>';
    if(deflectionText) h+='<h3>Angular Deflection</h3><p><strong>'+deflectionText+'</strong></p>';

    document.getElementById('detailedResults').innerHTML=h;
    document.getElementById('results').scrollIntoView({behavior:'smooth'});
}

function generateDiagrams(result,T){
    var x=result.x;var n=x.length;
    var torqueData=[];for(var i=0;i<n;i++)torqueData.push(T);

    document.getElementById('diagrams').innerHTML=
        '<h2>Horizontal Loading Diagram (HLD)</h2><div id="sfd_h"></div><div id="bmd_h"></div>'+
        '<h2>Vertical Loading Diagram (VLD)</h2><div id="sfd_v"></div><div id="bmd_v"></div>'+
        '<h2>Resultant Bending Moment</h2><div id="bmd_res"></div>'+
        '<h2>Torque Diagram</h2><div id="torque_d"></div>';

    var maxH=0,idxH=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_h[i])>maxH){maxH=Math.abs(result.BM_h[i]);idxH=i;}}
    var maxV=0,idxV=0;for(var i=0;i<n;i++){if(Math.abs(result.BM_v[i])>maxV){maxV=Math.abs(result.BM_v[i]);idxV=i;}}
    var maxR=0,idxR=0;for(var i=0;i<n;i++){if(result.BM_res[i]>maxR){maxR=result.BM_res[i];idxR=i;}}

    Plotly.newPlot('sfd_h',[{x:x,y:result.SF_h,type:'scatter',fill:'tozeroy',line:{color:'rgb(55,128,191)',width:3}}],
        {title:'Shear Force Diagram - Horizontal Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'SF (N)'}});

    Plotly.newPlot('bmd_h',[{x:x,y:result.BM_h,type:'scatter',fill:'tozeroy',line:{color:'rgb(219,64,82)',width:3}}],
        {title:'Bending Moment Diagram - Horizontal Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N-mm)'},
        annotations:[{x:x[idxH],y:result.BM_h[idxH],text:'M_max = '+result.BM_h[idxH].toFixed(0)+' N-mm',showarrow:true,arrowhead:2,ax:40,ay:-40}]});

    Plotly.newPlot('sfd_v',[{x:x,y:result.SF_v,type:'scatter',fill:'tozeroy',line:{color:'rgb(55,128,191)',width:3}}],
        {title:'Shear Force Diagram - Vertical Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'SF (N)'}});

    Plotly.newPlot('bmd_v',[{x:x,y:result.BM_v,type:'scatter',fill:'tozeroy',line:{color:'rgb(219,64,82)',width:3}}],
        {title:'Bending Moment Diagram - Vertical Plane',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N-mm)'},
        annotations:[{x:x[idxV],y:result.BM_v[idxV],text:'M_max = '+result.BM_v[idxV].toFixed(0)+' N-mm',showarrow:true,arrowhead:2,ax:40,ay:-40}]});

    Plotly.newPlot('bmd_res',[{x:x,y:result.BM_res,type:'scatter',fill:'tozeroy',line:{color:'rgb(148,0,211)',width:3}}],
        {title:'Resultant BM = sqrt(M_h^2 + M_v^2)',xaxis:{title:'Distance (mm)'},yaxis:{title:'BM (N-mm)'},
        annotations:[{x:x[idxR],y:maxR,text:'M_max = '+maxR.toFixed(0)+' N-mm',showarrow:true,arrowhead:2,ax:40,ay:-40}]});

    Plotly.newPlot('torque_d',[{x:x,y:torqueData,type:'scatter',fill:'tozeroy',line:{color:'rgb(50,171,96)',width:3}}],
        {title:'Torque Diagram',xaxis:{title:'Distance (mm)'},yaxis:{title:'Torque (N-mm)'}});
}
