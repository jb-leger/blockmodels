
dispatcher <- function(membership_name,membership_init,model_name,network){
	.Call( "dispatcher",
           membership_name,
           membership_init,
           model_name,
           network,
          PACKAGE = "lsbm" )
}

