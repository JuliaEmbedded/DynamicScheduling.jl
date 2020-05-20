#top level module where all user-facing functions are exported

module DynamicScheduling

################ exports ############################

################ using ##############################
using SSATools
using Core.Compiler
using InteractiveUtils
using LightGraphs
using Juno

################ type includes ######################
include("EC_types.jl")

################ types ##############################
#generate the elastic circuit - may move this to a separate file
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
    return ctr.val
end
function dec(ctr::counter)
    ctr.val-=1
    return ctr.val
end

function rst(ctr::counter)
    ctr.val=0
end

################ function includes ##################
#include("EC_DOT_printer.jl") #TODO workout the dependency/sub-module julia way of doing this

################ EC helper functions ################
function get_concrete_types(t::DataType, ret_arr::Array{DataType, 1})
        t_arr = InteractiveUtils.subtypes(t)

        if length(t_arr) == 0
                push!(ret_arr, t)
        else
                ret = vcat([get_concrete_types(type) for type in t_arr]...)
                append!(ret_arr, ret)
        end
        return ret_arr
end

get_concrete_types(t::DataType) = get_concrete_types(t, DataType[])

function convert_operator(node::SSATools.CDFGNode, node_idx::Int, inst_cntrs::Dict{DataType, counter}, arg_len::Int, node_len::Int, cnsts::Vector{ECconstant}) #TODO add to Base.convert
    predComps = Int[]
    input1Type=nothing
    input2Type=nothing
    int_types = get_concrete_types(Integer)
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
            if type ∉ int_types
                error("only integer constants supported")
            end
            cnst_idx =  node_len + arg_len + length(cnsts)+1
            push!(cnsts, ECconstant(:cst_, node.bb, inc(inst_cntrs[ECconstant]), val, type, Int[0], Int[node_idx]))
            push!(predComps, cnst_idx)
            inc(inst_cntrs[Any])
        end
    end

    if node.op.name == :mul_int
        return mul_int(:mul_, node.bb, inc(inst_cntrs[mul_int]), input1Type, input2Type, node.type, 0.000, 4, 1, predComps, copy(node.dataSuccs)), cnsts
    elseif node.op.name == :sub_int
        return sub_int(:sub_, node.bb, inc(inst_cntrs[sub_int]), input1Type, input2Type, node.type, 1.693, 0, 1, predComps, copy(node.dataSuccs)), cnsts
    elseif node.op.name == :add_int
        return add_int(:add_, node.bb, inc(inst_cntrs[sub_int]), input1Type, input2Type, node.type, 1.693, 0, 1, predComps, copy(node.dataSuccs)), cnsts
    elseif node.op.name == :slt_int
        return slt_int(:icmp_, node.bb, inc(inst_cntrs[slt_int]), input1Type, input2Type, node.type, 1.530, 0, 1, predComps, copy(node.dataSuccs)), cnsts
    elseif node.op.name == :sle_int
        return sle_int(:icmp_, node.bb, inc(inst_cntrs[sle_int]), input1Type, input2Type, node.type, 1.530, 0, 1, predComps, copy(node.dataSuccs)), cnsts
    elseif node.op.name == :(===)
        if input1Type ∈ int_types && input2Type ∈ int_types
            return eq_int(:icmp_, node.bb, inc(inst_cntrs[eq_int]), input1Type, input2Type, node.type, 1.530, 0, 1, predComps, copy(node.dataSuccs)), cnsts
        else
            error("Unsupported comparison of non integers")
        end
    elseif node.op.name == :checked_sdiv_int
        sdiv_int(:sdiv_, node.bb, inc(inst_cntrs[sdiv_int]), input1Type, input2Type, node.type, 0.966, 36, 1, predComps, copy(node.dataSuccs)), cnsts
    else
        error("Unsupported operator: ", string(node.op.name))
    end
end

function connect_mux_merge!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter}, MC_idxs::Array{Int, 1})
    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if isa(cmpt, mux)
            if MC_idxs[cmpt.bbID] == 0
                error("control merge missing")
            end
            if length(cmpt.predBBs) != length(ec.components[MC_idxs[cmpt.bbID]].predBBs)
                error("control merge does not have the same number of BBs")
            end

            new_mux_preds = []
            for (pred_idx, pred_bb) in zip(ec.components[MC_idxs[cmpt.bbID]].predComps, ec.components[MC_idxs[cmpt.bbID]].predBBs)
                #deconstruct component preds based on bbs
                bb_idx = findall(isequal(pred_bb), cmpt.predBBs)[1]
                push!(new_mux_preds, cmpt.predComps[bb_idx])
            end
            cmpt.predComps = new_mux_preds
            cmpt.predBBs = copy(ec.components[MC_idxs[cmpt.bbID]].predBBs)
            cmpt.predCtrl = MC_idxs[cmpt.bbID]

            push!(ec.components[MC_idxs[cmpt.bbID]].succCtrls, cmpt_idx)
        end
    end
    return ec, inst_cntrs
end

#adds the control infrastructure
function add_controls!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter})
    #looking at the ec c++ lib, they add the control information to bbnodes as well, from the dot graph printer this info is not used - redundant?
    #TODO check if this information is worth including in bbnodes - for now, not worth it since it requires type change of vector

    #TODO implement merge fork combo - assuming fork and other passes have been run

    #add end component - attaches to all return components, update in future
    end_inputTypes = DataType[]
    end_predComps = Int[]
    push!(ec.components, exit_ctrl(:end_, 0, inc(inst_cntrs[exit_ctrl]), end_inputTypes, end_predComps, [0])) #assuming succ doesnt matter atm
    inc(inst_cntrs[Any]) #update the overall instance counter
    end_ctrl_idx = length(ec.components)

    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if isa(cmpt, return_op)
            push!(cmpt.succComps, end_ctrl_idx)
            push!(ec.components[end_ctrl_idx].inputTypes, cmpt.output1Type)
            push!(ec.components[end_ctrl_idx].predComps, cmpt_idx)
        end
    end

    #add start and control phi connections - TODO add control routing for phi nodes
    ctrl_idxs = []
    MC_idxs = Int[]
    for (bb_num, bb) in enumerate(ec.bbnodes)
        ctrl_idx=0
        push!(MC_idxs, 0)
        if bb_num == 1
            #add the start components - should only be in bb1
            push!(ec.components, entry(:start_, 1, inc(inst_cntrs[entry]), true, Core.Any, Int[0], [])) #assuming succ doesnt matter atm
            ctrl_idx = length(ec.components)
        else
            #check the predecessors
            if length(bb["ctrlPreds"]) <= 0
                error("Middle bb has no predecessors")
            else
                dataPreds = Int[]
                inputTypes = DataType[]
                predBBs = Int[]

                ctrl_idx = length(ec.components)+1
                if length(bb["ctrlPreds"]) == 1
                    push!(ec.components, merge_data(:phiC_a, bb_num, inc(inst_cntrs[merge_data]), inputTypes, Core.Any, dataPreds, [], 0.166))
                else #more than two pred bbs -> create merge_ctrl
                    push!(ec.components, merge_ctrl(:phiMC_a, bb_num, inc(inst_cntrs[merge_ctrl]), inputTypes, Core.Any, dataPreds, [], [], 0.166, predBBs))
                    MC_idxs[end] = ctrl_idx
                end
                inc(inst_cntrs[Any]) #update the overall instance counter
                for predbb in bb["ctrlPreds"]
                    push!(inputTypes, Any) #control just passes tokens afaik
                    if length(ec.bbnodes[predbb]["ctrlSuccs"]) == 1 && !isa(ec.components[ec.bbnodes[predbb]["stmt_idx"][end]], branch) # Shouldn't be activated now
                        #generate new control branch and hook it up to drive any constants for the implicit fall through
                        println("GEN NEW BRANCH_C - bb: ", predbb)
                        cnst_idx = length(ec.components) + 1
                        branch_idx = length(ec.components) + 2
                        push!(ec.components, ECconstant(:brCst_, predbb, inc(inst_cntrs[ECconstant]), true, Bool, [0], [branch_idx]))
                        push!(ec.components, branch(:branchC_, predbb, inc(inst_cntrs[branch]), Any, Any, Any, [0, cnst_idx], ctrl_idx, 0)) #constant branch, must be true?
                        inc(inst_cntrs[Any]) #update the overall instance counter
                        inc(inst_cntrs[Any]) #update the overall instance counter
                        push!(ec.bbnodes[predbb]["stmt_idx"], branch_idx) # update bbnode collection with implicit branch

                        #get the ctrl branches that link to this phiC
                        push!(dataPreds, ec.bbnodes[predbb]["stmt_idx"][end])
                        push!(predBBs, predbb)
                    else
                        #get the ctrl branches that link to this phiC
                        push!(dataPreds, ec.bbnodes[predbb]["stmt_idx"][end]) #TODO make using "end" more robust
                        push!(predBBs, predbb)

                        #update those ctrl branches with the bT and bF successors
                        if (predbb+1) == bb_num #if it's true, it goes straight through
                            ec.components[ec.bbnodes[predbb]["stmt_idx"][end]].branchT = ctrl_idx
                        else
                            ec.components[ec.bbnodes[predbb]["stmt_idx"][end]].branchF = ctrl_idx
                        end
                    end
                end
            end
        end
        push!(ctrl_idxs, ctrl_idx)
    end
    #println(ctrl_idxs)
    #return ec, inst_cntrs
    for (bb, ctrl_idx) in zip(ec.bbnodes, ctrl_idxs)
        #don't think I can do this in one loop through
        #TODO hook up control flow to branches, sels, memory components, etc...
        if length(bb["ctrlSuccs"]) > 0
            #this node has successors which require control flow

            #assumptions
            #that gotoifnots will have created branches - WRONG, implicit fall throughs
            #only two children of a basicblock when gotoifnots are concerned
            #true - straight through (current bb+1)
            #false - the other bb in ctrlSuccs
            #branch statement exists at the end of the: bb["stmt_idx"][end]

            ec.components[bb["stmt_idx"][end]].predComps[1] = ctrl_idx #add phi to control branch preds
            ec.components[ctrl_idx].output1Type = Core.Any #TODO check this is valid when control signals are actually used
            push!(ec.components[ctrl_idx].succComps, bb["stmt_idx"][end]) #add ctrl branch to phi succs

            #hook up constants
            if isa(ec.components[ec.components[bb["stmt_idx"][end]].predComps[2]], ECconstant) #if branch selLine is a constant
                ec.components[ec.components[bb["stmt_idx"][end]].predComps[2]].predComps[1] = ctrl_idx
                push!(ec.components[ctrl_idx].succComps, ec.components[bb["stmt_idx"][end]].predComps[2])
            end

        else #no branching for this bb
            #llvm has implicit return - TODO check if julia tir has the same issue, I don't think so? but maybe if nothing returned
            push!(ec.components, sink(:sink_, 0, inc(inst_cntrs[sink]), Any, Int[ctrl_idx]))
            inc(inst_cntrs[Any])
            push!(ec.components[ctrl_idx].succComps, lastindex(ec.components)) #correct successor of the current control node
        end
    end

    ec, inst_cntrs = connect_mux_merge!(ec, inst_cntrs, MC_idxs)

    return ec, inst_cntrs
end

function add_phis!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter}, cfg)
    debug_flag=false
    bbnodes_cp = deepcopy(ec.bbnodes)

    #setting the predecessors of the phi nodes to placeholders - interfers with live in analysis
    plhdr_idxs = []
    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if isa(cmpt, mux)
            erase_nodes = Int[]
            add_nodes = Int[]
            for (pred, pred_bb) in zip(cmpt.predComps, cmpt.predBBs)
                push!(ec.components, placeholder(:phi_placeholder_, pred_bb, [pred], [cmpt_idx]))
                plhr_idx = length(ec.components)
                #preserves successor order
                succ_idx = findall(isequal(cmpt_idx), ec.components[pred].succComps)[1] #should only be one
                ec.components[pred].succComps[succ_idx] = plhr_idx
                push!(add_nodes, plhr_idx)
                push!(erase_nodes, pred)

                #pred_bb_def = SSATools.get_bb_num(cfg, pred)
                pred_bb_def = ec.components[pred].bbID #should be equivalent - necessary because pred may have been auto generated (not be in CDFG)
                #add to stmt_idx in bbs - defs
                pushfirst!(ec.bbnodes[pred_bb]["stmt_idx"], plhr_idx)
                #edit the targets of the preds - uses
                phi_idx = findall(isequal(cmpt_idx), ec.bbnodes[pred_bb_def]["stmt_tgts"][pred]["tgt_idx"])
                #change the og for the placeholders
                ec.bbnodes[pred_bb_def]["stmt_tgts"][pred]["tgt_idx"][phi_idx[1]] = plhr_idx
                ec.bbnodes[pred_bb_def]["stmt_tgts"][pred]["tgt_bb"][phi_idx[1]] = pred_bb

                push!(ec.bbnodes[pred_bb]["stmt_tgts"], plhr_idx=>Dict("tgt_bb"=>[cmpt.bbID], "tgt_idx"=>[cmpt_idx]))
            end
            filter!(rem-> rem ∉ erase_nodes, cmpt.predComps) #remove preds
            union!(cmpt.predComps, add_nodes)
            append!(plhdr_idxs, add_nodes)
        end
    end

    #live-ness analysis
    bb_defs = [copy(bb["stmt_idx"]) for bb in ec.bbnodes]
    bb_uses = [[] for _ in 1:length(ec.bbnodes)]
    for (bb_num, bb) in enumerate(ec.bbnodes)
        for stmt in bb["stmt_idx"]
            for (tgt_bb, tgt_idx) in zip(bb["stmt_tgts"][stmt]["tgt_bb"], bb["stmt_tgts"][stmt]["tgt_idx"])
                if tgt_bb != bb_num && !isa(ec.components[tgt_idx], mux) # ignore IR phis
                    #if stmt ∉ bb_uses[tgt_bb]
                    #    push!(bb_uses[tgt_bb], stmt)
                    #end
                    union!(bb_uses[tgt_bb], stmt)
                end
            end
        end
    end
    if debug_flag
        println("defs: ", bb_defs)
        println("uses: ", bb_uses)
    end

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

    if debug_flag
        println("live_ins: ", bb_live_ins)
    end

    phinodes = [[], []]
    for (bbl_num, bbl) in enumerate(bb_live_ins)
        for stmt in bbl
            []
            phi_a = merge_data(:phi_a, bbl_num, inc(inst_cntrs[merge_data]), [ec.components[stmt].output1Type], ec.components[stmt].output1Type, [stmt], [], 0.000)
            push!(phinodes[2], phi_a)
            push!(ec.components, phi_a) #works like you would expect, both ref the one instance
            inc(inst_cntrs[Any])

            phi_idx = length(ec.components)
            push!(ec.components[stmt].succComps, phi_idx)
            push!(phinodes[1], phi_idx)

            pushfirst!(ec.bbnodes[bbl_num]["stmt_idx"], phi_idx)
            pushfirst!(bbnodes_cp[bbl_num]["stmt_idx"], phi_idx-length(plhdr_idxs)) #placeholders added before phis
        end
    end

    if debug_flag
        println(length(phinodes[1]))
        for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
            println(phi_idx, " - ", phi)
        end
    end

    #hook up the phis where relevant - connections made within a block
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

    if debug_flag
        println("##DEBUG1##")
        for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
            println(phi_idx, " - ", phi)
        end
    end

    #phis to preds of phis - this forms chains to carry
    #for data though the different control flows
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        for bb_pred in ec.bbnodes[phi.bbID]["ctrlPreds"]
            phi_def_tgt = phi.predComps[1]
            for (phi_pred_idx, phi_pred) in zip(phinodes[1], phinodes[2])
                #phi pred is from a predecessor bb of the current phi, this isnt the og bb, they have the same og def stmt as pred
                if phi_pred.bbID == bb_pred && bb_pred != ec.components[phi_def_tgt].bbID && phi_def_tgt == phi_pred.predComps[1]
                    push!(phi_pred.succComps, phi_idx)

                    push!(phi.predComps, phi_pred_idx)
                    push!(phi.inputTypes, ec.components[phi_pred_idx].output1Type)
                end
            end
        end
    end
    if debug_flag
        println("##DEBUG2##")
        for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
            println(phi_idx, " - ", phi)
        end
    end

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

                    if isa(ec.components[succ], merge_data) || isa(ec.components[succ], mux) #This could be a very long list
                        del_idx = findall(isequal(cmpt_idx), ec.components[succ].predComps)
                        for idx in del_idx
                            deleteat!(ec.components[succ].predComps, idx)
                            deleteat!(ec.components[succ].inputTypes, idx)
                        end
                    else
                        filter!(rem->rem != cmpt_idx, ec.components[succ].predComps)
                    end
                else
                    #add the relevant phi nodes to the liveouts for branching later
                    if isa(cmpt, merge_data) && (ec.components[succ].bbID != cmpt.bbID || isa(ec.components[succ], merge_data))
                        union!(ec.bbnodes[cmpt.bbID]["liveouts"], [cmpt_idx])
                        union!(bbnodes_cp[cmpt.bbID]["liveouts"], [cmpt_idx-length(plhdr_idxs)])
                    end
                end
            end
        end
    end

    #exchange multi input phis for mux
    for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
        if length(ec.components[phi_idx].predComps) > 1
            ec.components[phi_idx] = mux(:phiM_a, phi.bbID, inc(inst_cntrs[mux]), phi.inputTypes, phi.output1Type, phi.predComps, phi.succComps, 0, 0.366, [ec.components[pred].bbID for pred in phi.predComps])
            dec(inst_cntrs[merge_data])
        end
    end

    if debug_flag
        println("##DEBUG3##")
        for (phi_idx, phi) in zip(phinodes[1], phinodes[2])
            println(phi_idx, " - ", phi)
        end
    end

    #merge through the placeholders
    for plhdr in reverse(plhdr_idxs)
        plhdr_bb = ec.components[plhdr].bbID
        for (pred, succ) in zip(ec.components[plhdr].predComps, ec.components[plhdr].succComps)
            pred_bb = ec.components[pred].bbID
            succ_bb = ec.components[succ].bbID
            #update the pred with the succ
            succ_idx = findall(isequal(plhdr), ec.components[pred].succComps)[1]
            ec.components[pred].succComps[succ_idx] = succ
            #update the succ with the pred
            pred_idx = findall(isequal(plhdr), ec.components[succ].predComps)[1]
            ec.components[succ].predComps[pred_idx] = pred
        end

        #remove the placeholders
        deleteat!(ec.components, plhdr)
        #adjust component indices
        for cmpt in ec.components
            if isa(cmpt, branch) # has predcomps, branchT, branchF
                cmpt.predComps = [(cmpt_idx >= plhdr ? cmpt_idx-1 : cmpt_idx) for cmpt_idx in cmpt.predComps]
                cmpt.branchT = (cmpt.branchT >= plhdr ? cmpt.branchT-1 : cmpt.branchT)
                cmpt.branchF = (cmpt.branchF >= plhdr ? cmpt.branchF-1 : cmpt.branchF)
            elseif isa(cmpt, source) # has succComps
                cmpt.succComps = [(cmpt_idx >= plhdr ? cmpt_idx-1 : cmpt_idx) for cmpt_idx in cmpt.succComps]
            elseif isa(cmpt, sink) # has predComps
                cmpt.predComps = [(cmpt_idx >= plhdr ? cmpt_idx-1 : cmpt_idx) for cmpt_idx in cmpt.predComps]
            else #other components have .predComps AND predSuccs
                cmpt.predComps = [(cmpt_idx >= plhdr ? cmpt_idx-1 : cmpt_idx) for cmpt_idx in cmpt.predComps]
                cmpt.succComps = [(cmpt_idx >= plhdr ? cmpt_idx-1 : cmpt_idx) for cmpt_idx in cmpt.succComps]
            end
        end
    end
    ec.bbnodes = bbnodes_cp #resets this to original(no placeholders) plus phi-based alterations
    return ec, inst_cntrs
end

function add_branches!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter})
    for (bb_num, bb) in enumerate(ec.bbnodes)
        goto_cmpt = ec.components[bb["stmt_idx"][end]]
        if length(bb["ctrlSuccs"]) <= 0 #has not bbs afterwards so can't branch from here (must have a return)
            #sanity check
            if !isa(goto_cmpt, return_op)
                error("BB leaf does not finish with a return - is this legal?")
            end
            continue
        else
            #cnst_flag = false
            if !isa(goto_cmpt, branch)
                #implicit fall through or something went very wrong - way to check
                #only one succ bb and curr bb in pred bb of succ
                if length(bb["ctrlSuccs"]) > 1
                    error("last component is not a branch: ", goto_cmpt)
                end

                push!(ec.components, ECconstant(:brCst_, bb_num, inc(inst_cntrs[ECconstant]), true, Core.Bool, Int[0], Int[]))
                sel_line_pred = length(ec.components)
                #push!(constants, ECconstant(:cst_, bb_num, inc(inst_cntrs[ECconstant]), true, Core.Bool, Int[0], Int[]))
                inc(inst_cntrs[Any])
                #cnst_flag = true

                #add the inevitable control branch required so that there aren't duplicate constants
                push!(ec.components, branch(:branchC_, bb_num, inc(inst_cntrs[branch]), Any, Any, Any, [0, sel_line_pred], 0 , 0))
                inc(inst_cntrs[Any])
                push!(bb["stmt_idx"], length(ec.components))
            else
                #hunt for switching input
                sel_line_pred = goto_cmpt.predComps[2] #second input is select line
                if sel_line_pred == 0
                    error("branch has no select line")
                end
            end

            for cmpt_idx in bb["liveouts"]
                #ASSUMPTION - need a branch per variable
                branch_idx_dec1 = length(ec.components)
                input1Type = (isa(ec.components[cmpt_idx], ECconstant) ? ec.components[cmpt_idx].type : ec.components[cmpt_idx].output1Type)
                branchT = 0
                branchF = 0
                for (bb_succ) in bb["ctrlSuccs"] #ASSUMPTION - bbs max out branching to two bbs
                    #look through the bb_succ phi nodes to see if the cmpt_idx(liveout) is the same as a pred
                    for tgt_stmt in ec.bbnodes[bb_succ]["stmt_idx"]
                        pred_list = findall(isequal(cmpt_idx), ec.components[tgt_stmt].predComps)
                        if length(pred_list) <= 0
                            #cmpt not targeting here
                        elseif length(pred_list) == 1
                            if !(isa(ec.components[tgt_stmt], merge_data) || isa(ec.components[tgt_stmt], mux))
                                error("should only be applying to phi nodes")
                            end

                            (bb_succ == bb_num+1 ? branchT = tgt_stmt : branchF = tgt_stmt) #select the correct branch for the tgt
                            ec.components[tgt_stmt].predComps[pred_list[1]] = branch_idx_dec1+1 #set tgt pred as branch

                            #remove tgt_stmt from succ of live out (cmpt_idx)
                            #println("succs: ", ec.components[cmpt_idx].succComps, "tgt: ", tgt_stmt)
                            filter!(rem->rem != tgt_stmt, ec.components[cmpt_idx].succComps)

                            #add the branch to the live out succs
                            union!(ec.components[cmpt_idx].succComps, [branch_idx_dec1+1])
                        else
                            error("multi targets for one component")
                        end
                    end
                end
                push!(ec.components, branch(:branch_, bb_num, inc(inst_cntrs[branch]), input1Type, input1Type, input1Type, Int[cmpt_idx, sel_line_pred], branchT, branchF))
                inc(inst_cntrs[Any])
                #update sel line succ
                union!(ec.components[sel_line_pred].succComps, [length(ec.components)])
            end
        end
    end
    return ec, inst_cntrs
end

#adds the forks
function add_forks!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter})
    no_fork_list = [branch, sink, exit_ctrl, fork] #component types with allowable multiple successors OR should not be forked from
    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if typeof(cmpt) ∉ no_fork_list && length(cmpt.succComps) > 1
            #println(cmpt)
            #start the forking
            fork_succs = []
            fork_idx = length(ec.components)+1
            for succ in cmpt.succComps[1:end] #WARNING not sure why this required, maybe because list is edited as its traversed?
                #check that its valid (comes from same bb)
                if cmpt.bbID != ec.components[succ].bbID && ec.components[succ].bbID != 0
                    error("between block linking should not occur at this point")
                end
                #add succs to fork succs
                push!(fork_succs, succ)
                #filter!(rem->rem != succ, cmpt.succComps) #remove succ from cmpt succ
                #filter!(rem->rem!=cmpt_idx, ec.components[succ].predComps)
                rpl_idx = findall(isequal(cmpt_idx), ec.components[succ].predComps)[1]
                ec.components[succ].predComps[rpl_idx] = fork_idx #replace cmpt from succ preds
            end

            if isa(cmpt, ECconstant) #need to expand this for other weird types
                fork_type = cmpt.type
            else
                fork_type = cmpt.output1Type
            end
            fork_name = :fork_
            if cmpt.name ∈ [:phiC_a, :phiMC_a]
                fork_name = :forkC_
            end
            push!(ec.components, fork(fork_name, cmpt.bbID, inc(inst_cntrs[fork]), fork_type, 0, fork_type, [cmpt_idx], fork_succs))
            inc(inst_cntrs[Any])

            #add fork to cmpt succ
            cmpt.succComps = [fork_idx] #replace succs
            #[union!(ec.components[succ].predComps, fork_idx) for succ in fork_succs] # add the fork as a pred to the og cmpt succs
            #println(cmpt)
        end
        #gen the mux fork
        if isa(cmpt, merge_ctrl) && length(cmpt.succCtrls) > 1
            bitWidth = Int64(ceil(log2(length(cmpt.predComps)))) #bitwidth of condition signal

            fork_succs = []
            fork_idx = length(ec.components)+1
            for succ in cmpt.succCtrls[1:end] #WARNING not sure why this required, maybe because list is edited as its traversed?
                #check that its valid (comes from same bb)
                if cmpt.bbID != ec.components[succ].bbID && ec.components[succ].bbID != 0
                    error("between block linking should not occur at this point")
                end
                #add succs to fork succs
                push!(fork_succs, succ)
                #filter!(rem->rem != succ, cmpt.succCtrls) #remove succ from cmpt succ
                #filter!(rem->rem!=cmpt_idx, ec.components[succ].predComps)
                #rpl_idx = findall(isequal(cmpt_idx), ec.components[succ].predComps)[1]
                ec.components[succ].predCtrl = fork_idx #replace cmpt from succ preds
            end
            fork_type = Any
            #add fork to cmpt succ
            cmpt.succCtrls = [fork_idx] #replace succs
            push!(ec.components, fork(:forkMC_, cmpt.bbID, inc(inst_cntrs[fork]), fork_type, bitWidth, fork_type, [cmpt_idx], fork_succs))
            inc(inst_cntrs[Any])
        end#
    end
    return ec, inst_cntrs
end

function add_sinks!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter})
    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if isa(cmpt, branch)
            sink_idx = 0
            empty_branch_cnt=0
            if ec.components[cmpt_idx].branchT == 0
                empty_branch_cnt+=1
                #add a sink to the branch
                push!(ec.components, sink(:sink_, cmpt.bbID, inc(inst_cntrs[sink]), cmpt.input1Type, [cmpt_idx]))
                inc(inst_cntrs[Any])
                sink_idx = length(ec.components)
                ec.components[cmpt_idx].branchT = sink_idx
            end
            if ec.components[cmpt_idx].branchF == 0
                empty_branch_cnt+=1
                push!(ec.components, sink(:sink_, cmpt.bbID, inc(inst_cntrs[sink]), cmpt.input1Type, [cmpt_idx]))
                inc(inst_cntrs[Any])
                sink_idx = length(ec.components)
                ec.components[cmpt_idx].branchF = sink_idx
            end

            if empty_branch_cnt >= 2
                error("branch has no successors")
            end
        end
    end
    return ec, inst_cntrs
end

function add_sources!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter})
    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if isa(cmpt, ECconstant)
            if cmpt.predComps[1] == 0 # don't always want to spawn tokens for all branches
                push!(ec.components, source(:source_, cmpt.bbID, inc(inst_cntrs[source]), cmpt.type, [cmpt_idx]))
                cmpt.predComps[1] = length(ec.components)
                inc(inst_cntrs[Any])
            end
        end
    end
    return ec, inst_cntrs
end

#buffers dont affect circuit functionality, might add more than necessary
#TODO check if merge_ctrl ctrl lines need buffering
function add_buffers_basic!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter})
    for (cmpt_idx, cmpt) in enumerate(ec.components)
        if typeof(cmpt) ∈ [merge_ctrl, merge_data, mux] #phi types
            if length(ec.bbnodes[cmpt.bbID]["ctrlPreds"]) > 1

                if cmpt.name ∈ [:phiC_a, :phiMC_a]
                    name = :buffC_
                else
                    name = :buff_
                end
                succComps = Int[]
                push!(ec.components, buffer(name, cmpt.bbID, inc(inst_cntrs[buffer]), cmpt.output1Type, [cmpt_idx], succComps))
                buff_idx = length(ec.components)

                #post fork pass so components should only have 1 succComp
                if length(cmpt.succComps) > 1
                    error("too many successors: ", cmpt)
                end

                tgt_cmpt = cmpt.succComps[1]
                tgt_pred_idx = findall(isequal(cmpt_idx), ec.components[tgt_cmpt].predComps)[1]
                push!(succComps, tgt_cmpt) #add cmpt succ to buffer
                ec.components[tgt_cmpt].predComps[tgt_pred_idx] = buff_idx #replace cmpt idx with buff idx in succComp
                cmpt.succComps[1] = buff_idx #replace succ idx with buff idx in cmpt
            end
        end
    end
    return ec, inst_cntrs
end

function add_buffers_cycles_only!(ec::ElasticCircuit, inst_cntrs::Dict{DataType, counter}, cfg::Core.Compiler.CFG)
    #dependent on Johnson cycle algorithm - limited in size of graphs
    dg = SSATools.cfg_to_lightgraph(cfg)
    cycles = LightGraphs.simplecycles(dg)

    if length(cycles) > 0 #graph has cycles so need for buffers
        loop_starts = Int[]
        for cycle in cycles
            start_bb = min(cycle...)
            union!(loop_starts, [start_bb])
        end

        for (cmpt_idx, cmpt) in enumerate(ec.components)
            if typeof(cmpt) ∈ [merge_ctrl, merge_data, mux] && cmpt.bbID ∈ loop_starts #phi types && head of a loop
                if cmpt.name ∈ [:phiC_a, :phiMC_a]
                    name = :buffC_
                else
                    name = :buff_
                end
                succComps = Int[]
                push!(ec.components, buffer(name, cmpt.bbID, inc(inst_cntrs[buffer]), cmpt.output1Type, [cmpt_idx], succComps))
                buff_idx = length(ec.components)

                #post fork pass so components should only have 1 succComp
                if length(cmpt.succComps) > 1
                    error("too many successors: ", cmpt)
                end

                tgt_cmpt = cmpt.succComps[1]
                tgt_pred_idx = findall(isequal(cmpt_idx), ec.components[tgt_cmpt].predComps)[1]
                push!(succComps, tgt_cmpt) #add cmpt succ to buffer
                ec.components[tgt_cmpt].predComps[tgt_pred_idx] = buff_idx #replace cmpt idx with buff idx in succComp
                cmpt.succComps[1] = buff_idx #replace succ idx with buff idx in cmpt
            end
        end
    end
    return ec, inst_cntrs
end

################ EC generation #######################

function ElasticCircuit(cdfg::SSATools.CDFG; add_buff=:cycle)::ElasticCircuit
    #bbnodes = copy(cdfg.cfg.blocks)
    #TODO going to have to implemnet something similar for args
    bbnodes = [Dict("stmt_idx"=>Int[idx for idx in block.stmts],
                    "stmt_tgts"=>Dict(idx=>Dict("tgt_bb"=>Int[], "tgt_idx"=>Int[]) for idx in block.stmts),
                    "ctrlPreds"=>Int[pred for pred in block.preds],
                    "ctrlSuccs"=>Int[succ for succ in block.succs],
                    "liveins"=>Int[],
                    "liveouts"=>Int[])
                    for block in cdfg.cfg.blocks]

    cmpts = AbstractElasticComponent[]

    inst_cntrs = Dict(ec_type=>counter() for ec_type in get_concrete_types(AbstractElasticComponent)) #instances of the AbstractElasticComponents
    push!(inst_cntrs, Core.Any=>counter()) #all instances
    arg_len = length(cdfg.args)
    node_len = length(cdfg.nodes)

    #convert the existing cdfgnodes to enodes, add them to the component list - sort out indexing/ssa numbering here
    #manual conversion for the time being - use type inference in future
    constants = ECconstant[]
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
                    union!(bbnodes[pred_bb]["liveouts"], [ec_idx])
                    #update the current bb live INS with the val value
                    union!(bbnodes[pred_bb]["liveins"], [ec_idx])
                else
                    #llvm phi live ins can originate from the same bb - TODO verify this is the case for IR
                    #set its status as live in and live out
                    #add val to live in and live out list
                end
            #else
            #    error("constants not supported yet")
            end
        end

        if isa(node.op, Core.GlobalRef) #it must be an operator
            op_tmp, constants = convert_operator(node, node_idx, inst_cntrs, arg_len, node_len, constants)
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
                    if type ∉ get_concrete_types(Integer)
                        error("only integer constants supported")
                    end
                    cnst_idx =  node_len + arg_len + length(constants)+1
                    push!(constants, ECconstant(:cst_, node.bb, inc(inst_cntrs[ECconstant]), node.dataPreds[1][1], node.dataPreds[2][1], Int[0], Int[node_idx]))
                    predComps[1] = cnst_idx
                    inc(inst_cntrs[Any])
                end
            else
                println("Returning nothing is valid")
            end
            op_tmp = return_op(:ret_, node.bb, inst_cntrs[return_op].val+1, input1Type, output1Type, 0.0, 0, 1, predComps, Int[]) #successor left as 0 (undefined) add exit later
        elseif node.op == :gotoifnot
            #add placeholder component to maintain index linking - TODO add a pass that removes these index links or ignore them in printing
            #op_tmp = placeholder(:gottoifnot_, node.bb, copy(node.dataPreds[1]), copy(node.dataSuccs) ) # this will lose constants or function arguments

            #add a control branch - might not work in all cases, we'll see
            if node.dataPreds[4][1] # literal bool
                if !isa(node.dataPreds[1][1], Core.SlotNumber)
                    if type ∉ get_concrete_types(Integer)
                        error("only integer constants supported")
                    end
                    cnst_idx =  node_len + arg_len + length(constants)+1
                    push!(constants, ECconstant(:brCst_, node.bb, inc(inst_cntrs[ECconstant]), node.dataPreds[1][1], node.dataPreds[2][1], Int[0], Int[node_idx]))
                    sel_line_pred = cnst_idx
                    inc(inst_cntrs[Any])
                else
                    sel_line_pred = (val.id-1 + node_len)
                end
            else
                sel_line_pred = node.dataPreds[1][1]
            end
            op_tmp = branch(:branchC_, node.bb, inc(inst_cntrs[branch]), Core.Any, Core.Any, Core.Any, Int[0, sel_line_pred], 0, 0)
        elseif node.op == :goto #handle like gotoifnot but with a constant value on the select line - set to true
            cnst_idx =  node_len + arg_len + length(constants)+1
            push!(constants, ECconstant(:brCst_, node.bb, inc(inst_cntrs[ECconstant]), true, Core.Bool, Int[0], Int[node_idx]))
            sel_line_pred = cnst_idx
            inc(inst_cntrs[Any])
            op_tmp = branch(:branchC_, node.bb, inc(inst_cntrs[branch]), Core.Any, Core.Any, Core.Any, Int[0, sel_line_pred], 0, 0)
        elseif node.op == :phi
            #change to a mux - sort control flow out in control stage
            predComps = Int[]
            inputTypes = DataType[]

            for (val, type, pos, lit_bool, bb_def) in zip(node.dataPreds[1], node.dataPreds[2], node.dataPreds[3], node.dataPreds[4], node.ctrlPreds)
                if !(lit_bool)
                    push!(predComps, val)
                elseif isa(val, Core.SlotNumber)
                    push!(predComps, (val.id-1 + node_len))
                else
                    if type ∉ get_concrete_types(Integer)
                        error("only integer constants supported")
                    end
                    cnst_idx =  node_len + arg_len + length(constants)+1
                    push!(constants, ECconstant(:cst_, bb_def, inc(inst_cntrs[ECconstant]), val, type, Int[0], Int[node_idx]))
                    push!(predComps, cnst_idx)
                    inc(inst_cntrs[Any])

                    #TODO verify
                    pushfirst!(bbnodes[bb_def]["stmt_idx"], cnst_idx)
                    push!(bbnodes[bb_def]["stmt_tgts"], cnst_idx=>Dict("tgt_bb"=>Int[node.bb], "tgt_idx"=>Int[node_idx]))
                    union!(bbnodes[bb_def]["liveouts"], [cnst_idx])
                    union!(bbnodes[node.bb]["liveins"], [cnst_idx])
                end
                push!(inputTypes, type)
            end
            op_tmp = mux(:phi_, node.bb, inc(inst_cntrs[mux]), inputTypes, node.type, predComps, copy(node.dataSuccs), 0, 0.366, copy(node.ctrlPreds))
        elseif node.op == :nth #treat instance of nothing as a place holder for an implicit branch link
            cnst_idx =  node_len + arg_len + length(constants)+1
            push!(constants, ECconstant(:brCst_, node.bb, inc(inst_cntrs[ECconstant]), true, Core.Bool, Int[0], Int[node_idx]))
            sel_line_pred = cnst_idx
            inc(inst_cntrs[Any])
            op_tmp = branch(:branchC_, node.bb, inc(inst_cntrs[branch]), Core.Any, Core.Any, Core.Any, Int[0, sel_line_pred], 0, 0)
        else
            error("Unsupported cdfg op type")
        end
        inc(inst_cntrs[Any])
        push!(cmpts, op_tmp) #add the new component to the list
    end

    #start by adding the args to the component list - these are converted to entry nodes (not sure about control)
    for arg in cdfg.args
        #how to handle fork situation at this stage
        #TODO assuming block1, maybe change this in cdfg if only used in a later block
        arg_tmp = entry(arg.name, 1, inc(inst_cntrs[entry]), false, arg.type, Int[0], copy(arg.dataSuccs))
        inc(inst_cntrs[Any])
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

    #add the constants to the end of the component list
    append!(cmpts, constants)

    #initial ec without control, branches, phis, forks, sources, sinks etc.
    ec = ElasticCircuit(bbnodes, cmpts)

    #run the ec passes in the right order, logic passes
    ec, inst_cntrs = add_phis!(ec, inst_cntrs, cdfg.cfg)
    ec, inst_cntrs = add_branches!(ec, inst_cntrs)
    ec, inst_cntrs = add_controls!(ec, inst_cntrs)
    #simple passes
    ec, inst_cntrs = add_forks!(ec, inst_cntrs)
    ec, inst_cntrs = add_sinks!(ec, inst_cntrs)
    ec, inst_cntrs = add_sources!(ec, inst_cntrs)

    if add_buff == :basic #basic buffer pass, over-eagerly generates buffers
        ec, inst_cntrs = add_buffers_basic!(ec, inst_cntrs)
    elseif add_buff ==:cycle
        ec, inst_cntrs = add_buffers_cycles_only!(ec, inst_cntrs, cdfg.cfg)
    end

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
    return string(cmpt.name, ": unsupported cmpt printer")
end

function printDOT_cmpt(cmpt::entry)
    typesize = (cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8)) #bit size
    dot_str = ""
    dot_str*= "\"$(cmpt.name)$(cmpt.control ? "0" : "")\" [type = \"Entry\", "
    dot_str*= (cmpt.control ? "control = \"true\", " : "")
    dot_str*= "bbID = $(cmpt.bbID), in = \"in1:$typesize\", out = \"out1:$typesize\"];"
    return dot_str
end

function printDOT_cmpt(cmpt::exit_ctrl)
    dot_str = ""
    dot_str*= "\"end_$(cmpt.instNum)\" [type = \"Exit\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    for (ip_t_num, ip_t) in enumerate(cmpt.inputTypes)
        dot_str*= "in$ip_t_num:$(ip_t.size*8) "
    end
    dot_str*= "\", out = \"out1:$(cmpt.inputTypes[1].size*8)\"];" #still think this op_t pointless, change my mind
    return dot_str
end

function printDOT_cmpt(cmpt::branch)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Branch\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8)) in2?:1\", "
    dot_str*= "out = \"out1+:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8)) "
    dot_str*= "out2-:$(cmpt.output2Type == Core.Any ? 0 : (cmpt.output2Type.size*8))\"];"
    return dot_str
end

function printDOT_cmpt(cmpt::fork)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Fork\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \"in1:$(cmpt.input1Type == Core.Any ? cmpt.bitWidth : (cmpt.input1Type == Core.Bool ? 1 : cmpt.input1Type.size*8))\", "
    dot_str*= "out = \""
    for (op_num, op) in enumerate(cmpt.succComps)
        dot_str*= "out$op_num:$(cmpt.output1Type == Core.Any ? cmpt.bitWidth : (cmpt.output1Type == Core.Bool ? 1 : cmpt.output1Type.size*8)) "
    end
    dot_str*= "\"];"
    return dot_str
end

function printDOT_cmpt(cmpt::merge_data)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Merge\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    for (ip_t_num, ip_t) in enumerate(cmpt.inputTypes)
        dot_str*= "in$ip_t_num:$(ip_t == Core.Any ? 0 : (ip_t.size*8)) "
    end
    dot_str*= "\", out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay)];"
    return dot_str
end

function printDOT_cmpt(cmpt::merge_ctrl) #"phiC_8" [type = "CntrlMerge", bbID= 4, in = "in1:0 in2:0 ", out = "out1:0 out2?:1", delay=0.166];
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"CntrlMerge\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \""
    for (ip_t_num, ip_t) in enumerate(cmpt.inputTypes)
        dot_str*= "in$ip_t_num:$(ip_t == Core.Any ? 0 : (ip_t.size*8)) "
    end
    dot_str*= "\", out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8)) "
    bitWidth = Int64(ceil(log2(length(cmpt.predComps))))
    dot_str*= "out2?:$(bitWidth)\", "
    dot_str*= "delay = $(cmpt.delay)];"
    return dot_str
end

function printDOT_cmpt(cmpt::mux) #"phi_n5" [type = "Mux", bbID= 4, in = "in1?:1 in2:32 in3:32 ", out = "out1:32", delay=0.366];
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Mux\", "
    bitWidth = Int64(ceil(log2(length(cmpt.predComps))))
    dot_str*= "bbID = $(cmpt.bbID), in = \"in1?:$bitWidth "
    for (ip_t_num, ip_t) in enumerate(cmpt.inputTypes)
        dot_str*= "in$(ip_t_num+1):$(ip_t == Core.Any ? 0 : (ip_t.size*8)) "
    end
    dot_str*= "\", out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay)];"
    return dot_str
end

function printDOT_cmpt(cmpt::sink)
    dot_str = ""
    dot_str*= "\"sink_$(cmpt.instNum)\" [type = \"Sink\", "
    dot_str*= "bbID = $(cmpt.bbID), in = \"in1:0\"];" #assuming 0 - may need to include some type info in sink struct
    return dot_str
end

function printDOT_cmpt(cmpt::source)
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Source\", "
    dot_str*= "bbID = $(cmpt.bbID), out = \"out1:$(cmpt.type == Core.Any ? 0 : (cmpt.type.size*8))\"];" #assuming 0 - may need to include some type info in sink struct
    return dot_str
end

function printDOT_cmpt(cmpt::ECconstant) #"cst_0" [type = "Constant", bbID= 3, in = "in1:32", out = "out1:32", value = "0x00000001"];
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Constant\", "
    dot_str*= "bbID = $(cmpt.bbID), "
    dot_str*="in = \"in1:$(cmpt.type == Core.Any ? 0 : (cmpt.type == Core.Bool ? 1 : cmpt.type.size*8))\", "
    dot_str*="out = \"out1:$(cmpt.type == Core.Any ? 0 : (cmpt.type == Core.Bool ? 1 : cmpt.type.size*8))\", "
    dot_str*="value = \"0x$(string(cmpt.value, base=16))\"];" #assuming 0 - may need to include some type info in sink struct
    return dot_str
end

function printDOT_cmpt(cmpt::buffer) #"buffI_0" [type = "Buffer", bbID= 3, in = "in1:32", out = "out1:32"];
    dot_str = ""
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" [type = \"Buffer\", "
    dot_str*= "bbID = $(cmpt.bbID), "
    dot_str*="in = \"in1:$(cmpt.type == Core.Any ? 0 : (cmpt.type == Core.Bool ? 1 : cmpt.type.size*8))\", "
    dot_str*="out = \"out1:$(cmpt.type == Core.Any ? 0 : (cmpt.type == Core.Bool ? 1 : cmpt.type.size*8))\"];"
    return dot_str
end

function printDOT_cmpt(cmpt::return_op)
    dot_str = ""
    dot_str*= "\"ret_$(cmpt.instNum)\" [type = \"Operator\", "
    dot_str*= "bbID = $(cmpt.bbID), op = \"ret_op\", in = \""
    dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type == Core.Bool ? 1 : cmpt.input1Type.size*8))\", "
    dot_str*= "out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type == Core.Bool ? 1 : cmpt.output1Type.size*8))\", "
    dot_str*= "delay = $(cmpt.delay), latency = $(cmpt.latency), II = $(cmpt.II)];"
    return dot_str
end

macro add_operator_printers(op_ts::Symbol...)
    for op_t in op_ts
        op = ""
        if op_t ∈ [:slt_int, :eq_int]
            op = "icmp_"
        end
        if op == [:sle_int]
            op *= "ult"
        else
            op = split(string(op_t), "_")[1]
        end

        @eval begin
            function printDOT_cmpt(cmpt::$(op_t))
                dot_str = ""
                dot_str*= "\"$(cmpt.name)$(cmpt.instNum)\" [type = \"Operator\", "
                dot_str*= "bbID = $(cmpt.bbID), op = \"$($op)_op\", in = \""
                dot_str*= "in1:$(cmpt.input1Type == Core.Any ? 0 : (cmpt.input1Type.size*8)) "
                dot_str*= "in2:$(cmpt.input2Type == Core.Any ? 0 : (cmpt.input2Type.size*8))\", "
                dot_str*= "out = \"out1:$(cmpt.output1Type == Core.Any ? 0 : (cmpt.output1Type == Core.Bool ? 1 : cmpt.output1Type.size*8))\", "
                dot_str*= "delay = $(cmpt.delay), latency = $(cmpt.latency), II = $(cmpt.II)];"
                return dot_str
            end
        end
    end
end

@add_operator_printers mul_int sub_int add_int sdiv_int sle_int slt_int eq_int

################ link printers #######################
#TODO fix the colours, these are dictated by the targets as well as shooters - separate helper func for ease
function printDOT_link(cmpt_idx::Int, cmpt::AbstractElasticComponent, cmpts::Vector{AbstractElasticComponent}, tab_num::Int) #template?
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= string(cmpt.name, ": unsupported link printer")
    return dot_str
end

function get_idx(tgt, list; only_one=true)
    if only_one #get the only one, fail if there are more
        idx_arr = findall(isequal(tgt), list)
        if length(idx_arr) > 1
            error("Connecting twice not supported")
        end
        return idx_arr[1]
    else
        return idx_arr = findall(isequal(tgt), list)
    end
end

function printDOT_link(cmpt_idx::Int, cmpt::entry, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(cmpt.name)$(cmpt.control ? "0" : "")\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "
    idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
    dot_str*= "[color = \"$(cmpt.control ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::exit_ctrl, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    #has no successors
end

function printDOT_link(cmpt_idx::Int, cmpt::branch, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    for (out_num, br) in enumerate([cmpt.branchT, cmpt.branchF]) #TODO - may need to sort out adding succ vec for branches, hopefully not
        dot_str = ""
        for n in 1:tab_num
            dot_str *= "\t"
        end
        dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
        dot_str*= "\"$(string(cmpts[br].name, cmpts[br].instNum))\" "

        idx = get_idx(cmpt_idx, cmpts[br].predComps)
        colour = "red"
        if cmpt.name == :branchC_
            colour = "gold3"
        elseif cmpt.name == :branch_
            colour = "blue"
        end

        dot_str*= "[color = \"$(colour)\", minlen = 3, from = \"out$(out_num)\", to = \"in$(isa(cmpts[br], mux) ? idx+1 : idx)\"];" #minlen seems constant
        return dot_str
    end
end

function printDOT_link(cmpt_idx::Int, cmpt::fork, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    for (out_num, succ) in enumerate(cmpt.succComps)
        dot_str = ""
        for n in 1:tab_num
            dot_str *= "\t"
        end
        dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
        dot_str*= "\"$(string(cmpts[succ].name, cmpts[succ].instNum))\" "

        #finds the refs to the current component in the target component
        idx_arr = get_idx(cmpt_idx, cmpts[succ].predComps, only_one=false)
        if length(idx_arr) > 1
            error("Connecting twice not supported")
        elseif length(idx_arr) < 1 && !isa(cmpts[succ], mux)
            error("not connected: ", cmpt)
        end
        colour = (cmpt.name == :forkC_ ? "gold3" : (isa(cmpts[succ], mux) ? "green" : "red"))
        dot_str*= "[color = \"$colour\", from = \"out$(out_num)\", to = \"in$(isa(cmpts[succ], mux) ? "1" : idx_arr[1])\"];"
        return dot_str
    end
end

function printDOT_link(cmpt_idx::Int, cmpt::merge_data, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    #finds the refs to the current component in the target component
    idx_arr = findall(isequal(cmpt_idx), cmpts[cmpt.succComps[1]].predComps)
    if length(idx_arr) > 1
        error("Connecting twice not supported")
    end

    dot_str*= "[color = \"$(cmpt.name == :phiC_a ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::merge_ctrl, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    #finds the refs to the current component in the target component
    idx_arr = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps, only_one=false)
    if length(idx_arr) > 1
        error("Connecting twice not supported")
    elseif length(idx_arr) < 1
        error("not connected")
    end

    dot_str*= "[color = \"$(cmpt.name == :phiMC_a ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx_arr[1])\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str

    #mux condition driver
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succCtrls[1]].name, cmpts[cmpt.succCtrls[1]].instNum))\" "
    dot_str*= "[color = \"green\", from = \"out2\", to = \"in1\"];$(length(cmpt.succCtrls) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::mux, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::sink, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    #has no successors
end

function printDOT_link(cmpt_idx::Int, cmpt::source, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::ECconstant, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
    dot_str*= "[color = \"$(cmpt.name == :brCst_ ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::buffer, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
    dot_str*= "[color = \"$(cmpt.name == :buffC_ ? "gold3" : "red")\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

function printDOT_link(cmpt_idx::Int, cmpt::return_op, cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
    dot_str = ""
    for n in 1:tab_num
        dot_str *= "\t"
    end
    dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
    dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

    idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
    dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
    return dot_str
end

macro add_operator_linkers(op_ts::Symbol...)
    for op_t in op_ts
        @eval begin
            function printDOT_link(cmpt_idx::Int, cmpt::$(op_t), cmpts::Vector{AbstractElasticComponent}, tab_num::Int)
                dot_str = ""
                for n in 1:tab_num
                    dot_str *= "\t"
                end
                dot_str*= "\"$(string(cmpt.name, cmpt.instNum))\" -> "
                dot_str*= "\"$(string(cmpts[cmpt.succComps[1]].name, cmpts[cmpt.succComps[1]].instNum))\" "

                idx = get_idx(cmpt_idx, cmpts[cmpt.succComps[1]].predComps)
                dot_str*= "[color = \"red\", from = \"out1\", to = \"in$(idx)\"];$(length(cmpt.succComps) > 1 ? "FOR FORKS SAKE" : "")"
                return dot_str
            end
        end
    end
end

@add_operator_linkers mul_int sub_int add_int sdiv_int sle_int slt_int eq_int

################### EC Printer ######################
#REPL print is component list (for debugging)
Base.show(io::IO, ::MIME"text/plain", ec::ElasticCircuit) = begin
    for cmpt in ec.components
        println(io, cmpt)
    end
end

#print() is the graph print
Base.show(io::IO, ec::ElasticCircuit) = begin
    #preamble
    println(io, "Digraph G {")
    println(io, "\tsplines=spline;")
    #println(io, "//DHLS version: 0.1.1\" [shape = \"none\" pos = \"20,20!\"]")

    bbnode_stmts = [ Int[] for i in ec.bbnodes] # containers to store the statements in a bbnode
    #print the components
    for (cmpt_idx,cmpt) in enumerate(ec.components)
        line = "\t\t" * printDOT_cmpt(cmpt)
        println(io, line)

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
        println(io, "\tsubgraph cluster_", string(bb_num-1), " {")
        println(io, "\tcolor = \"darkgreen\";")
        println(io, "\t\tlabel = \"block", string(bb_num), "\";")
        for cmpt_idx in bb
            line = printDOT_link(cmpt_idx, ec.components[cmpt_idx], ec.components, tab_num)
            if !isa(line, Nothing)
                println(io, line)
            end
        end
        println(io, "\t}")
    end
    #postamble
    println(io, "}")
end

#original interactive printing for Juno debugging
Juno.render(i::Juno.Inline, ec::ElasticCircuit) = Juno.render(i, Juno.defaultrepr(ec))

#################### DOT Flow #######################
function dot_from_f(@noinline(f), args, build_path)
        #process function
        f_tir = code_typed(f, args)[1]
        f_cdfg = SSATools.get_cdfg(f_tir.first)
        f_ec = DynamicScheduling.ElasticCircuit(f_cdfg)

        #write folder
        func_name = string(f)
        rm("build_path$(func_name)_sim/", force=true, recursive=true)
        mkpath(build_path * "build/$(func_name)_sim")

        graph_path = build_path * "build/$(func_name)_sim/$(func_name)_graph.dot"
        dot_f = open(graph_path, "w")
        print(dot_f, f_ec)
        close(dot_f)
        return graph_path
end

end #DynamicScheduling
