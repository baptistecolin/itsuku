with System.Storage_Elements;

package Phis is
   type Index is new Natural;
   subtype Phi_Seed is System.Storage_Elements.Storage_Array (1 .. 4);

   function Phi ( Seed : System.Storage_Elements.Storage_Array;
                  I : Index ) 
   return Index with 
      Pre  => Seed'Length = 4,
      Post => (0<=Phi'Result) and (Phi'Result<I);

end Phis;
