with Interfaces;

package Custom_Types is
   type Index is range 0 .. 2**63-1;
   for Index'Size use 64;
   type Byte is new Interfaces.Unsigned_8;
   type Byte_Array is array (Natural range <>) of Byte;
end Custom_Types;
