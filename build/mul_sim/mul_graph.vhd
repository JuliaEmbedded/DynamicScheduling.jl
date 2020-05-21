-- ==============================================================
-- Generated by Dot2Vhdl ver. 0.21
-- File created: Thu May 21 17:54:35 2020

-- ==============================================================
library IEEE; 
use IEEE.std_logic_1164.all; 
use IEEE.numeric_std.all; 
use work.customTypes.all; 
-- ==============================================================
entity mul_graph is 
port (
	clk:  in std_logic;
	rst:  in std_logic;
	start_in:  in std_logic_vector (0 downto 0);
	start_valid:  in std_logic;
	start_ready:  out std_logic;
	end_out:  out std_logic_vector (31 downto 0);
	end_valid:  out std_logic;
	end_ready:  in std_logic;
	a_din : in std_logic_vector (31 downto 0);
	a_valid_in : in std_logic;
	a_ready_out : out std_logic;
	b_din : in std_logic_vector (31 downto 0);
	b_valid_in : in std_logic;
	b_ready_out : out std_logic;
	c_din : in std_logic_vector (31 downto 0);
	c_valid_in : in std_logic;
	c_ready_out : out std_logic);
end;

architecture behavioral of mul_graph is 

	signal mul_0_clk : std_logic;
	signal mul_0_rst : std_logic;
	signal mul_0_dataInArray_0 : std_logic_vector(31 downto 0);
	signal mul_0_dataInArray_1 : std_logic_vector(31 downto 0);
	signal mul_0_pValidArray_0 : std_logic;
	signal mul_0_pValidArray_1 : std_logic;
	signal mul_0_readyArray_0 : std_logic;
	signal mul_0_readyArray_1 : std_logic;
	signal mul_0_nReadyArray_0 : std_logic;
	signal mul_0_validArray_0 : std_logic;
	signal mul_0_dataOutArray_0 : std_logic_vector(31 downto 0);

	signal mul_1_clk : std_logic;
	signal mul_1_rst : std_logic;
	signal mul_1_dataInArray_0 : std_logic_vector(31 downto 0);
	signal mul_1_dataInArray_1 : std_logic_vector(31 downto 0);
	signal mul_1_pValidArray_0 : std_logic;
	signal mul_1_pValidArray_1 : std_logic;
	signal mul_1_readyArray_0 : std_logic;
	signal mul_1_readyArray_1 : std_logic;
	signal mul_1_nReadyArray_0 : std_logic;
	signal mul_1_validArray_0 : std_logic;
	signal mul_1_dataOutArray_0 : std_logic_vector(31 downto 0);

	signal ret_0_clk : std_logic;
	signal ret_0_rst : std_logic;
	signal ret_0_dataInArray_0 : std_logic_vector(31 downto 0);
	signal ret_0_pValidArray_0 : std_logic;
	signal ret_0_readyArray_0 : std_logic;
	signal ret_0_nReadyArray_0 : std_logic;
	signal ret_0_validArray_0 : std_logic;
	signal ret_0_dataOutArray_0 : std_logic_vector(31 downto 0);

	signal a_clk : std_logic;
	signal a_rst : std_logic;
	signal a_dataInArray_0 : std_logic_vector(31 downto 0);
	signal a_pValidArray_0 : std_logic;
	signal a_readyArray_0 : std_logic;
	signal a_nReadyArray_0 : std_logic;
	signal a_validArray_0 : std_logic;
	signal a_dataOutArray_0 : std_logic_vector(31 downto 0);

	signal b_clk : std_logic;
	signal b_rst : std_logic;
	signal b_dataInArray_0 : std_logic_vector(31 downto 0);
	signal b_pValidArray_0 : std_logic;
	signal b_readyArray_0 : std_logic;
	signal b_nReadyArray_0 : std_logic;
	signal b_validArray_0 : std_logic;
	signal b_dataOutArray_0 : std_logic_vector(31 downto 0);

	signal c_clk : std_logic;
	signal c_rst : std_logic;
	signal c_dataInArray_0 : std_logic_vector(31 downto 0);
	signal c_pValidArray_0 : std_logic;
	signal c_readyArray_0 : std_logic;
	signal c_nReadyArray_0 : std_logic;
	signal c_validArray_0 : std_logic;
	signal c_dataOutArray_0 : std_logic_vector(31 downto 0);

	signal end_0_clk : std_logic;
	signal end_0_rst : std_logic;
	signal end_0_dataInArray_0 : std_logic_vector(31 downto 0);
	signal end_0_pValidArray_0 : std_logic;
	signal end_0_readyArray_0 : std_logic;
	signal end_0_nReadyArray_0 : std_logic;
	signal end_0_validArray_0 : std_logic;
	signal end_0_dataOutArray_0 : std_logic_vector(31 downto 0);
	signal end_0_validArray_1 :  std_logic;
	signal end_0_dataOutArray_1 :  std_logic_vector (31 downto 0);
	signal end_0_nReadyArray_1 :  std_logic;

	signal start_0_clk : std_logic;
	signal start_0_rst : std_logic;
	signal start_0_dataInArray_0 : std_logic_vector(0 downto 0);
	signal start_0_pValidArray_0 : std_logic;
	signal start_0_readyArray_0 : std_logic;
	signal start_0_nReadyArray_0 : std_logic;
	signal start_0_validArray_0 : std_logic;
	signal start_0_dataOutArray_0 : std_logic_vector(0 downto 0);

	signal sink_0_clk : std_logic;
	signal sink_0_rst : std_logic;
	signal sink_0_dataInArray_0 : std_logic_vector(0 downto 0);
	signal sink_0_pValidArray_0 : std_logic;
	signal sink_0_readyArray_0 : std_logic;

begin


	mul_0_clk <= clk;
	mul_0_rst <= rst;
	mul_1_pValidArray_0 <= mul_0_validArray_0;
	mul_0_nReadyArray_0 <= mul_1_readyArray_0;
	mul_1_dataInArray_0 <= mul_0_dataOutArray_0;

	mul_1_clk <= clk;
	mul_1_rst <= rst;
	ret_0_pValidArray_0 <= mul_1_validArray_0;
	mul_1_nReadyArray_0 <= ret_0_readyArray_0;
	ret_0_dataInArray_0 <= mul_1_dataOutArray_0;

	ret_0_clk <= clk;
	ret_0_rst <= rst;
	end_0_pValidArray_0 <= ret_0_validArray_0;
	ret_0_nReadyArray_0 <= end_0_readyArray_0;
	end_0_dataInArray_0 <= ret_0_dataOutArray_0;

	a_clk <= clk;
	a_rst <= rst;
	a_dataInArray_0 <= a_din;
	a_pValidArray_0 <= start_valid;
	mul_0_pValidArray_0 <= a_validArray_0;
	a_nReadyArray_0 <= mul_0_readyArray_0;
	mul_0_dataInArray_0 <= a_dataOutArray_0;

	b_clk <= clk;
	b_rst <= rst;
	b_dataInArray_0 <= b_din;
	b_pValidArray_0 <= start_valid;
	mul_0_pValidArray_1 <= b_validArray_0;
	b_nReadyArray_0 <= mul_0_readyArray_1;
	mul_0_dataInArray_1 <= b_dataOutArray_0;

	c_clk <= clk;
	c_rst <= rst;
	c_dataInArray_0 <= c_din;
	c_pValidArray_0 <= start_valid;
	mul_1_pValidArray_1 <= c_validArray_0;
	c_nReadyArray_0 <= mul_1_readyArray_1;
	mul_1_dataInArray_1 <= c_dataOutArray_0;

	end_0_clk <= clk;
	end_0_rst <= rst;
	end_valid <= end_0_validArray_0;
	end_out <= end_0_dataOutArray_0;
	end_0_nReadyArray_0 <= end_ready;

	start_0_clk <= clk;
	start_0_rst <= rst;
	start_0_pValidArray_0 <= start_valid;
	start_ready <= start_0_readyArray_0;
	sink_0_pValidArray_0 <= start_0_validArray_0;
	start_0_nReadyArray_0 <= sink_0_readyArray_0;
	sink_0_dataInArray_0 <= start_0_dataOutArray_0;

	sink_0_clk <= clk;
	sink_0_rst <= rst;

mul_0: entity work.mul_op(arch) generic map (2,1,32,32)
port map (
	clk => mul_0_clk,
	rst => mul_0_rst,
	dataInArray(0) => mul_0_dataInArray_0,
	dataInArray(1) => mul_0_dataInArray_1,
	pValidArray(0) => mul_0_pValidArray_0,
	pValidArray(1) => mul_0_pValidArray_1,
	readyArray(0) => mul_0_readyArray_0,
	readyArray(1) => mul_0_readyArray_1,
	nReadyArray(0) => mul_0_nReadyArray_0,
	validArray(0) => mul_0_validArray_0,
	dataOutArray(0) => mul_0_dataOutArray_0
);

mul_1: entity work.mul_op(arch) generic map (2,1,32,32)
port map (
	clk => mul_1_clk,
	rst => mul_1_rst,
	dataInArray(0) => mul_1_dataInArray_0,
	dataInArray(1) => mul_1_dataInArray_1,
	pValidArray(0) => mul_1_pValidArray_0,
	pValidArray(1) => mul_1_pValidArray_1,
	readyArray(0) => mul_1_readyArray_0,
	readyArray(1) => mul_1_readyArray_1,
	nReadyArray(0) => mul_1_nReadyArray_0,
	validArray(0) => mul_1_validArray_0,
	dataOutArray(0) => mul_1_dataOutArray_0
);

ret_0: entity work.ret_op(arch) generic map (1,1,32,32)
port map (
	clk => ret_0_clk,
	rst => ret_0_rst,
	dataInArray(0) => ret_0_dataInArray_0,
	pValidArray(0) => ret_0_pValidArray_0,
	readyArray(0) => ret_0_readyArray_0,
	nReadyArray(0) => ret_0_nReadyArray_0,
	validArray(0) => ret_0_validArray_0,
	dataOutArray(0) => ret_0_dataOutArray_0
);

a: entity work.start_node(arch) generic map (1,1,32,32)
port map (
	clk => a_clk,
	rst => a_rst,
	dataInArray(0) => a_dataInArray_0,
	pValidArray(0) => a_pValidArray_0,
	readyArray(0) => a_readyArray_0,
	nReadyArray(0) => a_nReadyArray_0,
	validArray(0) => a_validArray_0,
	dataOutArray(0) => a_dataOutArray_0
);

b: entity work.start_node(arch) generic map (1,1,32,32)
port map (
	clk => b_clk,
	rst => b_rst,
	dataInArray(0) => b_dataInArray_0,
	pValidArray(0) => b_pValidArray_0,
	readyArray(0) => b_readyArray_0,
	nReadyArray(0) => b_nReadyArray_0,
	validArray(0) => b_validArray_0,
	dataOutArray(0) => b_dataOutArray_0
);

c: entity work.start_node(arch) generic map (1,1,32,32)
port map (
	clk => c_clk,
	rst => c_rst,
	dataInArray(0) => c_dataInArray_0,
	pValidArray(0) => c_pValidArray_0,
	readyArray(0) => c_readyArray_0,
	nReadyArray(0) => c_nReadyArray_0,
	validArray(0) => c_validArray_0,
	dataOutArray(0) => c_dataOutArray_0
);

end_0: entity work.end_node(arch) generic map (1,0,1,32,32)
port map (
	clk => end_0_clk,
	rst => end_0_rst,
	dataInArray(0) => end_0_dataInArray_0,
	pValidArray(0) => end_0_pValidArray_0,
	readyArray(0) => end_0_readyArray_0,
	dataOutArray(0) => end_0_dataOutArray_0,
	validArray(0) => end_0_validArray_0,
	nReadyArray(0) => end_0_nReadyArray_0
);

start_0: entity work.start_node(arch) generic map (1,1,1,1)
port map (
	clk => start_0_clk,
	rst => start_0_rst,
	dataInArray(0) => start_0_dataInArray_0,
	pValidArray(0) => start_0_pValidArray_0,
	readyArray(0) => start_0_readyArray_0,
	nReadyArray(0) => start_0_nReadyArray_0,
	validArray(0) => start_0_validArray_0,
	dataOutArray(0) => start_0_dataOutArray_0
);

sink_0: entity work.sink(arch) generic map (1,0,1,32)
port map (
	clk => sink_0_clk,
	rst => sink_0_rst,
	dataInArray(0) => sink_0_dataInArray_0,
	pValidArray(0) => sink_0_pValidArray_0,
	readyArray(0) => sink_0_readyArray_0
);

end behavioral; 
