with System.Storage_Elements;
with Custom_Types; use Custom_Types;

package Phis is
   subtype Phi_Seed is Byte_Array (1 .. 4);

   function Phi ( Seed : Phi_Seed;
                  I : Index ) 
   return Index with 
      Pre  => Seed'Length = 4,
      Post => (0<=Phi'Result) and (Phi'Result<I);

private
   
   function Seed_To_Int ( Seed : Phi_Seed ) return Index
      with Post => Seed_To_Int'Result < 2**32;

end Phis;
