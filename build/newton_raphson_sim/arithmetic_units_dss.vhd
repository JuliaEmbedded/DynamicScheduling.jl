-----------------------------------------------------------------------
-- int div, version 0.0
-----------------------------------------------------------------------

Library IEEE;
use IEEE.std_logic_1164.all;

entity int_div is

    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        pvalid : IN STD_LOGIC;
        nready : IN STD_LOGIC;
        valid, ready : OUT STD_LOGIC;
        din0 : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
        din1 : IN STD_LOGIC_VECTOR(31 DOWNTO 0);
        dout : OUT STD_LOGIC_VECTOR(31 DOWNTO 0));
end entity;

architecture arch of int_div is

    component intop_sdiv_32ns_32ns_32_36_1_div is
    generic (
        in0_WIDTH   : INTEGER :=32;
        in1_WIDTH   : INTEGER :=32;
        out_WIDTH   : INTEGER :=32);
    port (
        clk         : in  STD_LOGIC;
        reset       : in  STD_LOGIC;
        ce          : in  STD_LOGIC;
        dividend    : in  STD_LOGIC_VECTOR(in0_WIDTH-1 downto 0);
        divisor     : in  STD_LOGIC_VECTOR(in1_WIDTH-1 downto 0);
        quot        : out STD_LOGIC_VECTOR(out_WIDTH-1 downto 0);
        remd        : out STD_LOGIC_VECTOR(out_WIDTH-1 downto 0));
    end component;

    signal q0, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13, q14, q15, q16, q17, q18, q19, q20, q21, q22, q23, q24, q25, q26, q27, q28, q29, q30, q31, q32, q33, q34: std_logic;
    signal remd : STD_LOGIC_VECTOR(31 downto 0);

begin
    intop_sdiv_32ns_32ns_32_36_1_div_U1 :  component intop_sdiv_32ns_32ns_32_36_1_div
    port map (
        clk => clk,
        reset => reset,
        ce => nready,
        dividend => din0,
        divisor => din1,
        quot => dout,
        remd => remd);

    ready<= nready;

       process (clk) is
       begin
          if rising_edge(clk) then  
            if (reset = '1') then
                q0 <= pvalid;
                q1 <= '0';
                q2 <= '0';
                q3 <= '0';
                q4 <= '0';
                q5 <= '0';
                q6 <= '0';
                q7 <= '0';
                q8 <= '0';
                q9 <= '0';
                q10 <= '0';
                q11 <= '0';
                q12 <= '0';
                q13 <= '0';
                q14 <= '0';
                q15 <= '0';
                q16 <= '0';
                q17 <= '0';
                q18 <= '0';
                q19 <= '0';
                q20 <= '0';
                q21 <= '0';
                q22 <= '0';
                q23 <= '0';
                q24 <= '0';
                q25 <= '0';
                q26 <= '0';
                q27 <= '0';
                q28 <= '0';
                q29 <= '0';
                q30 <= '0';
                q31 <= '0';
                q32 <= '0';
                q33 <= '0';
                q34 <= '0';
            elsif (nready='1') then
                q0 <= pvalid;
                q1 <= q0;
                q2 <= q1;
                q3 <= q2;
                q4 <= q3;
                q5 <= q4;
                q6 <= q5;
                q7 <= q6;
                q8 <= q7;
                q9 <= q8;
                q10 <= q9;
                q11 <= q10;
                q12 <= q11;
                q13 <= q12;
                q14 <= q13;
                q15 <= q14;
                q16 <= q15;
                q17 <= q16;
                q18 <= q17;
                q19 <= q18;
                q20 <= q19;
                q21 <= q20;
                q22 <= q21;
                q23 <= q22;
                q24 <= q23;
                q25 <= q24;
                q26 <= q25;
                q27 <= q26;
                q28 <= q27;
                q29 <= q28;
                q30 <= q29;
                q31 <= q30;
                q32 <= q31;
                q33 <= q32;
                q34 <= q33;
             end if;
          end if;
       end process;

       valid <= q34;

end architecture;

Library IEEE;
use IEEE.std_logic_1164.all;
use ieee.numeric_std.all;
use work.customTypes.all;

entity div_op is
Generic (
 INPUTS: integer; OUTPUTS: integer; DATA_SIZE_IN: integer; DATA_SIZE_OUT: integer
);
    port (
        clk : IN STD_LOGIC;
        rst : IN STD_LOGIC;
        pValidArray : IN std_logic_vector(1 downto 0);
        nReadyArray : in std_logic_vector(0 downto 0);
        validArray : out std_logic_vector(0 downto 0);
        readyArray : OUT std_logic_vector(1 downto 0);
        dataInArray : in data_array (1 downto 0)(DATA_SIZE_IN-1 downto 0); 
        dataOutArray : out data_array (0 downto 0)(DATA_SIZE_OUT-1 downto 0));
end entity;

architecture arch of div_op is

signal alu_out: STD_LOGIC_VECTOR (DATA_SIZE_OUT-1 downto 0);
signal alu_valid, alu_ready : STD_LOGIC;

signal join_valid : STD_LOGIC;
signal join_ReadyArray : std_logic_vector(1 downto 0);
signal nready_tmp : std_logic;
begin 

    join_write_temp:   entity work.join(arch) generic map(2)
            port map( pValidArray,  --pValidArray
                      alu_ready,     --nready                    
                      join_valid,         --valid          
                      readyArray);   --readyarray 

     divide: entity work.int_div (arch)
            port map (clk, rst,
                     join_valid,     --pvalid,
                     nready_tmp,         --nready,
                     alu_valid, --valid,
                     alu_ready, --ready,
                     dataInArray(0),           --din0
                     dataInArray(1),           --din1
                     alu_out);  --dout

dataOutArray(0) <= alu_out;
validArray(0) <= alu_valid;
nready_tmp <= nReadyArray(0);

end architecture;
