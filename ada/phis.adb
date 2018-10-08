with System.Storage_Elements;
with Custom_Types; use Custom_Types;

package body Phis is
   function Phi ( Seed : Phi_Seed;
                  I : Index)
   return Index is
      R : Index := I-1;
   begin
      return I-1; -- TODO : write actual body
   end;

end Phis;
