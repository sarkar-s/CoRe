"""
Contains commands for Reactome database query using Cypher
"""

command_set = {}
command_set['events'] = "MATCH(n:TopLevelPathway {displayName:'#',speciesName:'Homo sapiens'})-[:hasEvent*]->(rle:ReactionLikeEvent) RETURN DISTINCT rle"
command_set['sub events'] = "MATCH(n:Pathway {stId:'#',speciesName:'Homo sapiens'})-[:hasEvent*]->(rle:ReactionLikeEvent) RETURN rle"

command_set['input'] = {}
command_set['output'] = {}

command_set['input']['coarse'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:input*]->(pe:PhysicalEntity) RETURN DISTINCT pe"
command_set['output']['coarse'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:output*]->(pe:PhysicalEntity) RETURN DISTINCT pe"

command_set['input']['medium'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:input|hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe"
command_set['output']['medium'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:output|hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe"

#command_set['input']['fine'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:input|hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe"
#command_set['output']['fine'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:output|hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe"

command_set['input']['fine'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:input|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe"
command_set['output']['fine'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:output|hasComponent|hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe"

command_set['superpathways'] = "MATCH(r:ReactionLikeEvent {stId:'#'})<-[:hasEvent*]-(sp:Pathway) RETURN sp.displayName"
command_set['regulation'] = "MATCH(r:ReactionLikeEvent {stId:'#'})-[:regulatedBy]->(re:Regulation) RETURN DISTINCT re"
command_set['regulator'] = "MATCH(n:# {stId:'%'})-[:regulator]->(pe:PhysicalEntity) RETURN DISTINCT pe"
command_set['reverse'] = "MATCH(r:Reaction {stId:'#'})-[:reverseReaction*]->(re:Reaction) RETURN re"
command_set['entity'] = "MATCH(ge:GenomeEncodedEntity {stId:'#'}) RETURN ge"
command_set['reference_gene'] = "MATCH(ge:EntityWithAccessionedSequence{stId:'#'})-[:referenceEntity]->(re:ReferenceEntity) RETURN re"
command_set['subpathways'] = "MATCH(tp:TopLevelPathway{displayName:'#',speciesName:'Homo sapiens'})-[:hasEvent]->(p:Pathway) RETURN DISTINCT p.displayName"

command_set['members'] = "MATCH(p:PhysicalEntity {stId:'#',speciesName:'Homo sapiens'})-[:hasMember|hasCandidate*]->(pe:PhysicalEntity) RETURN DISTINCT pe.stId"
