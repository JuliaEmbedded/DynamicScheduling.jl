#top level module where all user-facing functions are exported

module DynamicScheduling

################ exports ############################

################ using ##############################
using SSATools
using Core.Compiler

################ type includes ######################
include("EC_types.jl")

################ types ##############################
#generate the elastic circuit - may move this to a separate file

################ function includes ##################
#include("EC_DOT_printer.jl") #TODO workout the dependency/sub-module julia way of doing this

################ EC helper functions ################
function convert_operator(node::SSATools.CDFGNode, inst_cnt::Int, arg_len::Int, node_len::Int) #TODO add to Base.convert
    predComps = Int[]
    input1Type=nothing
    input2Type=nothing
    for (val, type, pos, lit_bool) in zip(node.dataPreds[1], node.dataPreds[2], node.dataPreds[3], node.dataPreds[4])
        if pos == 1
            input1Type = type
        elseif pos == 2
            input2Type = type
        else
            error("position fail - more than two args")
        end

        if !(lit_bool)
            push!(predComps, val)
        elseif isa(val, Core.SlotNumber)
            push!(predComps, (val.id-1 + node_len))
        else
            #TODO handle constant values - best way might be to assume const are after args and return some sort of flag?
            error("unsupported constant values")
        end
    end

    if node.op.name == :mul_int
        return mul_int(:mul_, node.bb, inst_cnt+1, input1Type, input2Type, node.type, 0.000, 4, 1, predComps, copy(node.dataSuccs))
    elseif node.op.name == :sub_int
        return sub_int(:sub_, node.bb, inst_cnt+1, input1Type, input2Type, node.type, 1.693, 0, 1, predComps, copy(node.dataSuccs))
    elseif node.op.name == :slt_int
        return slt_int(:icmp_, node.bb, inst_cnt+1, input1Type, input2Type, node.type, 1.530, 0, 1, predComps, copy(node.dataSuccs))
    else
        error("Unsupported operator: ", string(node.op.name))
    end
end

#adds the control infrastructure
function add_controls!(ec::ElasticCircuit, inst_cnt::Int)
    #looking at the ec c++ lib, they add the control information to bbnodes as well, from the dot graph printer this info is not used - redundant?
    #TODO check if this information is worth including in bbnodes - for now, not worth it since it requires type change of vector

    #TODO implement merge fork combo - assuming fork and other passes have been run

    #add end component - attaches to all return components, update in future
    end_inputTypes = DataType[]
    end_predComps = Int[]
    push!(ec.components, exit_ctrl(:end_, 0, inst_cnt+1, end_inputTypes, end_predComps, [0])) #assuming succ doesnt matter atm
    inst_cnt+=1 #update the overall instance counter
    end_ctrl_idx = length(ec.components)

    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if isa(cmpt, return_op)
            push!(cmpt.succComps, end_ctrl_idx)
            push!(ec.components[end_ctrl_idx].inputTypes, cmpt.output1Type)
            push!(ec.components[end_ctrl_idx].predComps, cmpt_idx)
        end
    end

    #add start and control phi connections - TODO add control routing for phi nodes
    for (bb_num, bb) in enumerate(ec.bbnodes)
        if bb_num == 1
            #add the start components - should only be in bb1
            push!(ec.components, entry(:start_, 1, inst_cnt+1, true, Core.Any, Int[0], [])) #assuming succ doesnt matter atm
        else
            #check the predecessors
            if length(bb["ctrlPreds"]) <= 0
                error("Middle bb has no predecessors")
            else
                dataPreds = Int[]
                inputTypes = DataType[]
                for predbb in bb["ctrlPreds"]
                    #get the ctrl branches that link to this phiC
                    push!(inputTypes, Any) #control just passes tokens afaik
                    push!(dataPreds, ec.bbnodes[predbb]["stmt_idx"][end]) #TODO make using "end" more robust

                    #update those ctrl branches with the bT and bF successors
                    if (predbb+1) == bb_num #if it's true, it goes straight through
                        ec.components[ec.bbnodes[predbb]["stmt_idx"][end]].branchT = length(ec.components)+1
                    else
                        ec.components[ec.bbnodes[predbb]["stmt_idx"][end]].branchF = length(ec.components)+1
                    end
                end
            end
            push!(ec.components, merge_ctrl(:phiC_, bb_num, inst_cnt+1, inputTypes, Core.Any, dataPreds, [], 0.166))
        end

        inst_cnt+=1 #update the overall instance counter
        ctrl_idx = length(ec.components)

        #TODO hook up control flow to branches, sels, memory components, etc...
        if length(bb["ctrlSuccs"]) > 0
            #this node has successors which require control flow
            #error("branching not implemented yet")

            #assumptions
            #that gotoifnots will have created branches - TODO work out how gotos should be handled
            #only two children of a basicblock when gotoifnots are concerned
            #true - straight through (current bb+1)
            #false - the other bb in ctrlSuccs
            #branch statement exists at the end of the: bb["stmt_idx"][end]

            #push!(ec.components[bb["stmt_idx"][end]].predComps, ctrl_idx) #branchC
            ec.components[bb["stmt_idx"][end]].predComps[1] = ctrl_idx
            ec.components[ctrl_idx].output1Type = Core.Any #TODO check this is valid when control signals are actually used
            push!(ec.components[ctrl_idx].succComps, bb["stmt_idx"][end]) #lastindex(ec.components)

        else #no branching for this bb
            #llvm has implicit return - TODO check if julia tir has the same issue, I don't think so? but maybe if nothing returned
            push!(ec.components, sink(:sink_, 0, inst_cnt+1, Int[ctrl_idx]))
            inst_cnt+=1
            push!(ec.components[ctrl_idx].succComps, lastindex(ec.components)) #correct successor of the current control node
        end

        #phi node creation here?

    end
end

#adds the forks
function add_forks!(ec::ElasticCircuit, inst_cnt::Int)
end

function add_sinks!(ec::ElasticCircuit, inst_cnt::Int)
end

function add_phis!(ec::ElasticCircuit, inst_cnt::Int)
    #TODO handle the mess that is IR phi nodes
    #live-ness analysis
    bb_defs = [copy(bb["stmt_idx"]) for bb in ec.bbnodes]
    bb_uses = [[] for _ in 1:length(ec.bbnodes)]
    for (bb_num, bb) in enumerate(ec.bbnodes)
        for stmt in bb["stmt_idx"]
            for (tgt_bb, tgt_idx) in zip(bb["stmt_tgts"][stmt]["tgt_bb"], bb["stmt_tgts"][stmt]["tgt_idx"])
                if tgt_bb != bb_num
                    if stmt ∉ bb_uses[tgt_bb]
                        push!(bb_uses[tgt_bb], stmt)
                    end
                end
            end
        end
    end
    #println("defs: ", bb_defs)
    #println("uses: ", bb_uses)

    bb_live_ins = [[] for _ in 1:length(ec.bbnodes)]

    mod_flag = true
    while mod_flag
        mod_flag = false
        for (bb_num, bb) in enumerate(ec.bbnodes)
            liveouts = Int[]
            [Base.union!(liveouts, bb_live_ins[bb_succ]) for bb_succ in bb["ctrlSuccs"]]

            diff = Base.setdiff(liveouts, bb_defs[bb_num])
            tmp = Base.union(bb_uses[bb_num], diff)

            if length(tmp) > length(bb_live_ins[bb_num])
                bb_live_ins[bb_num] = tmp
                mod_flag = true
            end
        end
    end

    #println("live_ins: ", bb_live_ins)

    phinodes = [[], []]
    for (bbl_num, bbl) in enumerate(bb_live_ins)
        for stmt in bbl
            phi_n = merge_ctrl(:phi_n, bbl_num, inst_cnt+1, [Any], Any, [stmt], [], 0.000)
            push!(phinodes[2], phi_n)
            push!(ec.components, phi_n) #works like you would expect, both ref the one instance
            inst_cnt+=1

            phi_idx = length(ec.components)
            push!(ec.components[stmt].succComps, phi_idx)
            push!(phinodes[1], phi_idx)
        end
    end
    #=println(length(phinodes[1]))
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        println(phi_idx, " - ", phi)
    end=#

    #hook up the phis where relevant - connections made withing a block
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        phi_def_idx = phi.predComps[1] #origin of the data phi node is transferring
        for (tgt_idx, tgt) in enumerate(ec.components)
            if (tgt.bbID == phi.bbID) && (tgt.bbID != ec.components[phi_def_idx].bbID) && (tgt != phi) && phi_def_idx ∈ tgt.predComps
                filter!(rem->rem != tgt_idx, ec.components[phi_def_idx].succComps) #remove the tgt/cmpt from succs of original component def
                filter!(rem->rem != phi_def_idx, tgt.predComps) #remove the original component def from preds of tgt/cmpt
                #put the phi nodes in between the connections
                push!(tgt.predComps, phi_idx) #add the phi as a pred to tgt/cmpt
                push!(phi.succComps, tgt_idx) #add the tgt/cmpt as a succ of phi
            end
        end
    end

    #=println("##DEBUG1##")
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        println(phi_idx, " - ", phi)
    end=#

    #phis to preds of phis - this forms chains to carry data though the different control flows
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        for bb_pred in ec.bbnodes[phi.bbID]["ctrlPreds"]
            phi_def_tgt = phi.predComps[1]
            for (phi_pred_idx, phi_pred) in zip(phinodes[1], phinodes[2])
                #phi pred is from a predecessor bb of the current phi, this isnt the og bb, they have the same og def stmt as pred
                if phi_pred.bbID == bb_pred && bb_pred != ec.components[phi_def_tgt].bbID && phi_def_tgt == phi_pred.predComps[1]
                    push!(phi_pred.succComps, phi_idx)
                    push!(phi.predComps, phi_pred_idx)
                end
            end
        end
    end
    #=println("##DEBUG2##")
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        println(phi_idx, " - ", phi)
    end=#

    for (cmpt_idx, cmpt) in enumerate(ec.components)
        #println(cmpt)
        if isa(cmpt, branch)
            #ignore
        else
            for succ in cmpt.succComps[1:end] #WARNING I have no idea why I need to be explicit about this array
                #println("cmpt_idx: ", cmpt_idx, " succ: ", succ, " succ_bb: ", ec.components[succ].bbID)
                if ec.components[succ].bbID ∉ ec.bbnodes[cmpt.bbID]["ctrlSuccs"] && cmpt.bbID != ec.components[succ].bbID
                    #yeet that out
                    filter!(rem->rem != succ, cmpt.succComps)
                    filter!(rem->rem != cmpt_idx, ec.components[succ].predComps)
                end
            end
        end
    end

    #=println("##DEBUG3##")
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        println(phi_idx, " - ", phi)
    end=#
end

mutable struct counter
    val::Int64

    function counter(val::Int)
        new_ctr = new()
        new_ctr.val = val
        return new_ctr
    end
end
counter() = counter(0)
function inc(ctr::counter)
    ctr.val+=1
end

function rst(ctr::counter)
    ctr.val=0
end

################ EC generation #######################

function ElasticCircuit(cdfg::SSATools.CDFG)::ElasticCircuit
    #bbnodes = copy(cdfg.cfg.blocks)
    #TODO going to have to implemnet something similar for args
    bbnodes = [Dict("stmt_idx"=>Int[idx for idx in block.stmts],
                    "stmt_tgts"=>Dict(idx=>Dict("tgt_bb"=>Int[], "tgt_idx"=>Int[]) for idx in block.stmts),
                    #"tgt_bb"=>Dict(idx=>Int[] for idx in block.stmts),
                    #"tgt_idx"=>Dict(idx=>Int[] for idx in block.stmts),
                    "ctrlPreds"=>Int[pred for pred in block.preds],
                    "ctrlSuccs"=>Int[succ for succ in block.succs],
                    "liveins"=>Int[],
                    "liveouts"=>Int[])
                    for block in cdfg.cfg.blocks]



    cmpts = AbstractElasticComponent[]

    inst_cnt = 0 #counts the component instances, TODO change to per component types in a dictionary or smth
    arg_len = length(cdfg.args)
    node_len = length(cdfg.nodes)

    #convert the existing cdfgnodes to enodes, add them to the component list - sort out indexing/ssa numbering here
    #manual conversion for the time being - use type inference in future
    for (node_idx, node) in enumerate(cdfg.nodes) #TODO try and come up with a more efficient way of doing this

        #only really care about pred information - I don't think it will be more efficient from successor direction
        for succ in node.dataSuccs
            push!(bbnodes[node.bb]["stmt_tgts"][node_idx]["tgt_idx"], copy(succ))
            push!(bbnodes[node.bb]["stmt_tgts"][node_idx]["tgt_bb"], SSATools.get_bb_num(cdfg.cfg, succ)) #find the appropriate bb
        end

        #compute bb live ins/outs
        for (val, type, pos, lit_bool) in zip(node.dataPreds[1], node.dataPreds[2], node.dataPreds[3], node.dataPreds[4])
            if !lit_bool || isa(val, Core.SlotNumber)#predecessor ssa val
                pred_bb = (isa(val, Core.SlotNumber) ? 1 : SSATools.get_bb_num(cdfg.cfg, val)) #if not int this bb
                ec_idx = (isa(val, Core.SlotNumber) ? val.id-1 + node_len : val)
                if pred_bb != node.bb
                    #update the pred bb live OUTS with the val idx
                    if ec_idx ∉ bbnodes[pred_bb]["liveouts"]
                        push!(bbnodes[pred_bb]["liveouts"], ec_idx)
                    end
                    #update the current bb live INS with the val value
                    if ec_idx ∉ bbnodes[node.bb]["liveins"]
                        push!(bbnodes[node.bb]["liveins"], ec_idx)
                    end
                else
                    #llvm phi live ins can originate from the same bb - TODO verify this is the case for IR
                    #set its status as live in and live out
                    #add val to live in and live out list
                end
            else
                error("constants not supported yet")
            end
        end

        if isa(node.op, Core.GlobalRef) #it must be an operator
            op_tmp = convert_operator(node, inst_cnt, arg_len, node_len)
        elseif node.op == :return #special kind of operator that doesnt use global ref
            input1Type=nothing
            output1Type=nothing
            predComps= Int[0]

            if length(node.dataPreds[1]) > 0
                input1Type = node.dataPreds[2][1]
                output1Type = node.dataPreds[2][1]

                if !node.dataPreds[4][1] #literal bool
                    predComps[1] = node.dataPreds[1][1]
                elseif isa(node.dataPreds[1][1], Core.SlotNumber)
                    predComps[1] = node.dataPreds[1][1].id-1 + node_len
                else
                    #handle constant values
                    error("unsupported constant values")
                end
            else
                println("Returning nothing is valid")
            end
            op_tmp = return_op(:ret_, node.bb, inst_cnt+1, input1Type, output1Type, 0.0, 0, 1, predComps, Int[]) #successor left as 0 (undefined) add exit later
        elseif node.op == :gotoifnot
            #add placeholder component to maintain index linking - TODO add a pass that removes these index links or ignore them in printing
            #op_tmp = placeholder(:gottoifnot_, node.bb, copy(node.dataPreds[1]), copy(node.dataSuccs) ) # this will lose constants or function arguments

            #add a control branch - might not work in all cases, we'll see
            if node.dataPreds[4][1]
                if !isa(node.dataPreds[1][1], Core.SlotNumber)
                    error("Constants not yet supported")
                else
                    sel_line_pred = (val.id-1 + node_len)
                end
            else
                sel_line_pred = node.dataPreds[1][1]
            end
            #op_tmp = branch(:branchC_, node.bb, inst_cnt+1, Core.Any, Core.Any, Core.Any, sel_line_pred, Int[], 0, 0)
            op_tmp = branch(:branchC_, node.bb, inst_cnt+1, Core.Any, Core.Any, Core.Any, Int[0, sel_line_pred], 0, 0)
        else
            error("Unsupported cdfg op type")
        end
        inst_cnt+=1
        push!(cmpts, op_tmp) #add the new component to the list
    end

    #start by adding the args to the component list - these are converted to entry nodes (not sure about control)
    for arg in cdfg.args
        #how to handle fork situation at this stage
        #TODO assuming block1, maybe change this in cdfg if only used in a later block
        arg_tmp = entry(arg.name, 1, inst_cnt+1, false, arg.type, Int[0], copy(arg.dataSuccs))
        inst_cnt+=1
        push!(cmpts, arg_tmp)

        arg_idx = length(cmpts)
        pushfirst!(bbnodes[1]["stmt_idx"], arg_idx) #arg idx ref
        push!(bbnodes[1]["stmt_tgts"], arg_idx=>Dict("tgt_bb"=>Int[], "tgt_idx"=>Int[]))

        for succ in arg.dataSuccs
            succ_bb = SSATools.get_bb_num(cdfg.cfg, succ)
            push!(bbnodes[1]["stmt_tgts"][arg_idx]["tgt_idx"], copy(succ))
            push!(bbnodes[1]["stmt_tgts"][arg_idx]["tgt_bb"], succ_bb) #find the appropriate bb
        end

        if arg_idx ∉ bbnodes[1]["liveins"]
            #println("ACTUALLY USEFUL")
            push!(bbnodes[1]["liveins"], arg_idx) #args are live ins of bb 1 ALWAYS (according to EC)
        end
    end

    #initial ec without control, branches, phis, forks, sources, sinks etc.
    ec = ElasticCircuit(bbnodes, cmpts)

    #run the ec passes in the right order
    add_phis!(ec, inst_cnt)
    add_controls!(ec, inst_cnt)

    #add_phis!(ec, inst_cnt)
    #add_forks!(ec, inst_cnt)
    #add_sinks!(ec, inst_cnt)

    #TODO some final validation before returning the ec
    return ec
end

#allows you to tap into the different levels of generation
ElasticCircuit(ci::Core.Compiler.CodeInfo) = (ci.inferred ? ElasticCircuit(SSATools.get_cdfg(ci)) : error("Lowered IR not supported at the moment"))
ElasticCircuit(ci_pair::Pair) = ElasticCircuit(ci_pair.first)
#input the function and the args of the function-> Tuple{Int64, etc}
ElasticCircuit(func, args::Tuple) = ElasticCircuit(code_typed(func, args)[1])

################ DOT printers ########################
#TODO change from println() to file storage
################ component printers ##################
function printDOT_cmpt(cmpt::AbstractElasticComponent) #template?
    println(cmpt.name, ": unsupported cmpt printer")
end

function printDOT_cmpt(cmpt::entry)
    typesize = (cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8)) #bit size
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Entry\", "
    dot_str*= (cmpt.control ? "control = \"true\", " : "")
    dot_str*= "bbID = $(cmpt.bbID), in = \"in1:$typesize\", out = \"out1:$typesize\"];"
    println(dot_str)
end

function printDOT_cmpt(cmpt::exit_ctrl)
    dot_str = ""
    dot_str*= "\"end_$(cmpt.instNum)\" [type = \"Exit\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    for (ip_t_num, ip_t) in enumerate(cmpt.inputTypes)
        dot_str*= "in$ip_t_num:$(ip_t.size*8) "
    end
    dot_str*= "\", out = \"out1:$(cmpt.inputTypes[1].size*8)\"];" #still think this op_t pointless, change my mind
    println(dot_str)
end

function printDOT_cmpt(cmpt::branch)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Branch\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8)) in2?:1\", "
    dot_str*= "out = \"out1+:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8)) "
    dot_str*= "out2-:$(cmpt.output2Type == Core.Any ? 0 : (cmpt.output2Type.size*8))\"];"
    println(dot_str)
end

function printDOT_cmpt(cmpt::merge_ctrl)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Merge\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    for (ip_t_num, ip_t) in enumerate(cmpt.inputTypes)
        dot_str*= "in$ip_t_num:$(ip_t == Core.Any ? 0 : (ip_t.size*8)) "
    end
    dot_str*= "\", out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay)];"
    println(dot_str)
end

function printDOT_cmpt(cmpt::sink)
    dot_str = ""
    dot_str*= "\"sink_$(cmpt.instNum)\" [type = \"Sink\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \"in1:0\"];" #assuming 0 - may need to include some type info in sink struct
    println(dot_str)
end

function printDOT_cmpt(cmpt::return_op)
    dot_str = ""
    dot_str*= "\"ret_$(cmpt.instNum)\" [type = \"Operator\", "
    dot_str*= "bbID = $(cmpt.bbID), op = \"ret_op\", in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8))\", "
    dot_str*= "out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay), latency = $(cmpt.latency), II = $(cmpt.II)];"
    println(dot_str)
end

function printDOT_cmpt(cmpt::mul_int)
    dot_str = ""
    dot_str*= "\"mul_$(cmpt.instNum)\" [type = \"Operator\", "
    dot_str*= "bbID = $(cmpt.bbID), op = \"mul_op\", in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8)) "
    dot_str*= "in2:$(cmpt.input2Type == Core.Any ? 0 : (cmpt.input2Type.size*8))\", "
    dot_str*= "out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay), latency = $(cmpt.latency), II = $(cmpt.II)];"
    println(dot_str)
end

function printDOT_cmpt(cmpt::sub_int)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Operator\", "
    dot_str*= "bbID = $(cmpt.bbID), op = \"sub_op\", in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8)) "
    dot_str*= "in2:$(cmpt.input2Type == Core.Any ? 0 : (cmpt.input2Type.size*8))\", "
    dot_str*= "out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay), latency = $(cmpt.latency), II = $(cmpt.II)];"
    println(dot_str)
end

function printDOT_cmpt(cmpt::slt_int)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Operator\", "
    dot_str*= "bbID = $(cmpt.bbID), op = \"icmp_slt_op\", in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8)) "
    dot_str*= "in2:$(cmpt.input2Type == Core.Any ? 0 : (cmpt.input2Type.size*8))\", "
    dot_str*= "out = \"out1:1\"" #$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", " #might need to be forced to 1 bit
    dot_str*= "delay = $(cmpt.delay), latency = $(cmpt.latency), II = $(cmpt.II)];"
    println(dot_str)
end

################ link printers #######################
#TODO fix the colours, these are dictated by the targets as well as shooters - separate helper func for ease
function printDOT_link(cmpt_idx::Int, cmpt::AbstractElasticComponent, cmpts::Vector{AbstractElasticComponent}, tab_num::Int) #template?
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    println(dot_str, cmpt.name, ": unsupported link printer")
end

function printDOT_link(cmpt_idx::Int, cmpt::entry, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps) #finds the refs to the current component in the target component
    if length(idx_arr) > 1
        error("Connecting twice not currently supported")
    end

    dot_str*= "[color = \"$(cmpt.control ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    println(dot_str)
end

function printDOT_link(cmpt_idx::Int, cmpt::exit_ctrl, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    #has no successors
end

function printDOT_link(cmpt_idx::Int, cmpt::branch, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    for br in [cmpt.branchT, cmpt.branchF] #TODO - may need to sort out adding succ vec for branches, hopefully not
        dot_str = ""
        for n in 1:tab_num
            dot_str *= "\t"
        end
        dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
        dot_str*= "\"$(string(cmpts[br].name, cmpts[br].instNum))\" "

        idx_arr = findall(isequal(cmpt_idx), cmpts[br].predComps) #finds the refs to the current component in the target component
        if length(idx_arr) > 1
            error("Connecting twice not currently supported")
        end

        dot_str*= "[color = \"$(cmpt.name == :branchC_ ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx_arr[1])\"];"
        println(dot_str)
    end
end

function printDOT_link(cmpt_idx::Int, cmpt::merge_ctrl, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps) #finds the refs to the current component in the target component
    if length(idx_arr) > 1
        error("Connecting twice not currently supported")
    end

    dot_str*= "[color = \"$(cmpt.name == :phiC_ ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    println(dot_str)
end

function printDOT_link(cmpt_idx::Int, cmpt::sink, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    #has no successors
end

function printDOT_link(cmpt_idx::Int, cmpt::return_op, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps) #finds the refs to the current component in the target component
    if length(idx_arr) > 1
        error("Connecting twice not currently supported")
    end

    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    println(dot_str)
end

function printDOT_link(cmpt_idx::Int, cmpt::mul_int, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps) #finds the refs to the current component in the target component
    if length(idx_arr) > 1
        error("Connecting twice not currently supported")
    end

    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    println(dot_str)
end

function printDOT_link(cmpt_idx::Int, cmpt::sub_int, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps) #finds the refs to the current component in the target component
    if length(idx_arr) > 1
        error("Connecting twice not currently supported")
    end

    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    println(dot_str)
end

function printDOT_link(cmpt_idx::Int, cmpt::slt_int, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps) #finds the refs to the current component in the target component
    if length(idx_arr) > 1
        error("Connecting twice not currently supported")
    end

    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    println(dot_str)
end

################ EC DOT Printer ######################
function printDOT(ec::ElasticCircuit)
    #preamble
    println("Digraph G {")
    println("\tsplines=spline;")
    println("//DHLS version: 0.1.1\" [shape = \"none\" pos = \"20,20!\"]")

    bbnode_stmts = [ Int[] for i in ec.bbnodes] # containers to store the statements in a bbnode
    #print the components
    for (cmpt_idx,cmpt) in enumerate(ec.components)
        print("\t\t")
        printDOT_cmpt(cmpt)

        if cmpt.bbID > 0
            push!(bbnode_stmts[cmpt.bbID], cmpt_idx)
        else
            #might need to change this bb assumption for blockless components
            push!(bbnode_stmts[1], cmpt_idx)
        end
    end
    #print the directed connections
    tab_num=2
    for (bb_num,bb) in enumerate(bbnode_stmts)
        println("\tsubgraph cluster_", string(bb_num-1), " {")
        println("\tcolor = \"darkgreen\";")
        println("\t\tlabel = \"block", string(bb_num), "\";")
        for cmpt_idx in bb
            printDOT_link(cmpt_idx, ec.components[cmpt_idx], ec.components, tab_num)
        end
        println("\t}")
    end

    #postamble
    println("}")
end

end #DynamicScheduling
