bool correctPhi         //在fvSolution文件夹的pimple子字典定义
(
    pimple.dict().getOrDefault("correctPhi", true)
);
 
bool checkMeshCourantNo
(
    pimple.dict().getOrDefault("checkMeshCourantNo", false)
);
 
bool moveMeshOuterCorrectors
(
    pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
);
 
 
bool massFluxInterpolation
(
    pimple.dict().getOrDefault("massFluxInterpolation", false)
);
 
bool adjustFringe
(
    pimple.dict().getOrDefault("oversetAdjustPhi", false)
);
 
bool ddtCorr
(
    pimple.dict().getOrDefault("ddtCorr", true)
);
