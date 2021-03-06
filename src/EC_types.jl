#top
abstract type AbstractElasticComponent end
#ec abstracts
abstract type AbstractOperator <: AbstractElasticComponent end
abstract type AbstractControl <: AbstractElasticComponent end
abstract type AbstractConnector <: AbstractElasticComponent end
#operator abstracts
abstract type AbstractMultiplier <: AbstractOperator end
abstract type AbstractDivider <: AbstractOperator end
abstract type AbstractAdder <: AbstractOperator end
abstract type AbstractSub <: AbstractOperator end
abstract type AbstractCompare <: AbstractOperator end

#concrete types - we're using "components" not nodes for elastic stuff, let's be consistent
mutable struct ElasticCircuit
    #bbnodes::Vector{Core.Compiler.BasicBlock} #BasicBlock members - stmts, preds, succs
    bbnodes::Vector{Any}
    components::Vector{AbstractElasticComponent} #data links are inside
end

#placeholders - component will always be replaced or removed
mutable struct placeholder <: AbstractElasticComponent
    #name::String #just take the type
    name::Symbol
    bbID::Int

    predComps::Vector{Int} # the idx of the phi pred its replacing
    succComps::Vector{Int} # phi mux idx
end

#=
struct gtin_placeholder <: AbstractElasticComponent
    name::Symbol # gotoifnot_
    bbID::Int
    predComps::Vector{Int} #driving control boolean (should only be one value afaik)
    succComps::Vector{Int} #shoud be zero ??

    #only potentially branched in bb1 - can't be sure when an arg will be used
    #so safest to assume it branches from bb1 - unsure how tokens spawned for arg control
    args::Vector{Int}

    bbstmts::Vector{Int}
end
=#

#elastic components
mutable struct fork <: AbstractConnector #maybe make mutable?
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    bitWidth::Int #ONLY used for linking merge_ctrl and phi mux nodes - zero otherwise
    output1Type::DataType #more than one output but they'll all be the same type

    predComps::Vector{Int} #input components (array positions) - should only be one
    succComps::Vector{Int} #output components (array positions)
end

mutable struct merge_data <: AbstractConnector
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    inputTypes::Vector{DataType}
    output1Type::DataType
    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)

    delay::Float32
end

mutable struct merge_ctrl <: AbstractConnector
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    inputTypes::Vector{DataType}
    output1Type::DataType
    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)

    succCtrls::Vector{Int}

    delay::Float32
    #helpful information
    predBBs::Vector{Int}
end

mutable struct mux <: AbstractConnector
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    inputTypes::Vector{DataType}
    output1Type::DataType #multiple types output currently not supported - type inference should help
    predComps::Vector{Int} #input components (array positions) first position is selecting input to route to op
    succComps::Vector{Int} #output component (array position)

    predCtrl::Int # should only be one

    delay::Float32 #0.366, might be just for two input

    #helpful information
    predBBs::Vector{Int}
end

mutable struct entry <: AbstractControl
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    #argNum::Int #order that arg appears in function entry
    control::Bool #false for args, true for control flow (not data)

    output1Type::DataType
    predComps::Vector{Int} #don't know why this is in original EC
    succComps::Vector{Int} #output component (array position)
end

mutable struct exit_ctrl <: AbstractControl
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    inputTypes::Vector{DataType}
    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #not sure why this is here in OGEC - zero is nc
end

mutable struct branch <: AbstractControl
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    output1Type::DataType #TODO fix - bit redundant to store the types as it just passes through
    output2Type::DataType

    #selLine::Int #the predecessor component that controls the branch - mutable as this may change
    predComps::Vector{Int} #input component (array value) - select line will be the second input
    branchT::Int #connection if selLine true
    branchF::Int #connection if selLine true
end


mutable struct source <: AbstractElasticComponent
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    type::DataType
    succComps::Vector{Int} #output component (array position)
end

mutable struct sink <: AbstractElasticComponent
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    type::DataType
    predComps::Vector{Int} #input components (array positions)
end

mutable struct ECconstant <: AbstractElasticComponent
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    value::Any
    type::DataType
    predComps::Vector{Int} #predecessors necessary for sauce
    succComps::Vector{Int} #output component (array position)
end

mutable struct buffer <: AbstractElasticComponent
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    type::DataType
    predComps::Vector{Int} #predecessors necessary for sauce
    succComps::Vector{Int} #output component (array position)
end

#operators
mutable struct return_op <: AbstractOperator #SISO for now
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    output1Type::DataType

    delay::Float32 #timed value
    latency::Int #cycles
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position) - zero if undefined
end

#mathematical operators
mutable struct mul_int <: AbstractMultiplier
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value
    latency::Int #cycles
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end

mutable struct sub_int <: AbstractSub
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value
    latency::Int #cycles
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end

mutable struct slt_int <: AbstractCompare # strictly less than, (first arg slt second arg)
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value
    latency::Int #cycles
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end

mutable struct eq_int <: AbstractCompare # equal - icmp_eq_op
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value - 1.530
    latency::Int #cycles - 0
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end

mutable struct sle_int <: AbstractCompare # strictly less than or equal, (first arg slt second arg)
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value
    latency::Int #cycles
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end

mutable struct add_int <: AbstractAdder #add_op
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value - 1.693
    latency::Int #cycles - 0
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end

mutable struct sdiv_int <: AbstractDivider
    name::Symbol
    bbID::Int #basic block number
    instNum::Int #the number of components of this type (prevents naming conflicts)

    input1Type::DataType
    input2Type::DataType
    output1Type::DataType

    delay::Float32 #timed value - 0
    latency::Int #cycles - 36
    II::Int #from what I've seen - always 1

    predComps::Vector{Int} #input components (array positions)
    succComps::Vector{Int} #output component (array position)
end
