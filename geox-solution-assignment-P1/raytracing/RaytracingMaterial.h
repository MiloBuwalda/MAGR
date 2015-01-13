//---------------------------------------------------------------------------
#ifndef RaytracingMaterialH
#define RaytracingMaterialH
//---------------------------------------------------------------------------
#include "PTypes.h"
//---------------------------------------------------------------------------



///
class RaytracingMaterial : public Persistent {
	GEOX_CLASS(RaytracingMaterial)
 public:
	Vector3f diffuseColor;
	Vector3f specularColor;
	Vector3f ambientColor;
	float32 reflectivity;
	float32 specularExponent;

	RaytracingMaterial();
};





#endif