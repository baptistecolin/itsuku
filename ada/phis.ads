with System.Storage_Elements;
with Custom_Types; use Custom_Types;

package Phis is
   type Index is new Natural;
   subtype Phi_Seed is System.Storage_Elements.Storage_Array (1 .. 4);

   function Phi ( Seed : Phi_Seed;
                  I : Index ) 
   return Index with 
      Pre  => Seed'Length = 4,
      Post => (0<=Phi'Result) and (Phi'Result<I);

end Phis;
