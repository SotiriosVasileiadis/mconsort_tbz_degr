
function Expander(oEl, oNotifier) {
    if (!oEl) oEl = document.body;
    var aElExpandHandler = document.querySelectorAll(".expand-handler");
    if (!aElExpandHandler) return;

    var sClassName = "e-hidden";
    for (var i = 0; i < aElExpandHandler.length; i++) {
        // prevent multiple 'addEvent' assignment
        if (utils.hasClass(aElExpandHandler[i], "_inited")) continue;
        utils.addEvent(aElExpandHandler[i], "click", x_Click);
        utils.addClass(aElExpandHandler[i], "_inited");
        var oElParent = x_FindParent(utils.getParent(aElExpandHandler[i]));
        if (utils.hasClass(oElParent, sClassName)) {
            utils.addClass(oElParent, sClassName);
        } else {
            utils.removeClass(oElParent, sClassName);
        }
    }

    function x_Click(e) {
        utils.preventDefault(e);
        var oElParent = utils.getParent(this);
        if (utils.hasClass(oElParent, "e-hidden")) {
            utils.addClass(oElParent, "showed");
            utils.removeClass(oElParent, "e-hidden");
        } else {
            utils.addClass(oElParent, "e-hidden");
            utils.removeClass(oElParent, "showed");
        }
        return false;
    }
 
    function x_FindParent(oEl) {
        while (oEl && !utils.hasClass(oEl, "expand")) {
            oEl = utils.getParent(oEl);
        }
        return oEl;
    }
}


//$Id: esummary.1.js 188644 2010-04-13 15:51:35Z sponomar $
/*
File: esummary.1.js
Class: Esummary
Send request to esummary and convert its XML output to javascript object.
Only <Item> element will be used.

<Item Name="Organism_Name" Type="String">Bacillus subtilis subsp. natto BEST195</Item>
will be converted to {..., Organism_Name: "Bacillus%20subtilis%20subsp.%20natto%20BEST195", ...}
*/


/*
Function: Esummary()
*/
function Esummary() {
    this.oCache = {};
}

/*
 ******************************************************************************************************
*/
Esummary.prototype.Lookup = function(oA) {
    var oThis = this;
    var sGet = this.Get(oA);
    
    var oProvider = new RemoteDataProvider("/entrez/eutils/esummary.fcgi?");
    oProvider.bAsync  = false;
    oProvider.onError = function(oObj) {}
    oProvider.onSuccess = function(oObj) {
        var n;
        try { n = oObj.responseXML.getElementsByTagName("Item"); } catch (e) { return; }
        var obj = {};
        for (var i = 0; i < n.length; i++) {
            var s = n[i].textContent;
            if (s == undefined) s = n[i].text;
            if (s > "") obj[n[i].getAttribute("Name")] = escape(s);
        }
        oThis.oCache[sGet] = obj;
        oThis.Process(oA, obj);
    }
    
    if (this.oCache[sGet]) {
        this.Process(oA, this.oCache[sGet]);
    } else {
        oProvider.Get(sGet);
        
    }
}


/*
 ******************************************************************************************************
*/
Esummary.prototype.Get = function(oA) {
    return "db=genomeprj&id=" + oA.getAttribute("rel");
}

/*
 ******************************************************************************************************
 this placeholder must be overwritten
*/
Esummary.prototype.Process = function(oA, obj) {
//    console.info(obj);
}


function Taxonomy() {};
Taxonomy.prototype  = new Esummary();
Taxonomy.prototype.Get = function(oA) {
    return "db=taxonomy&id=" + oA.getAttribute("rel");
}
Taxonomy.prototype.Process = function (oA, obj) {
    oA.innerHTML = unescape(obj.ScientificName);
}



/*
 *******************************************************************************
 User found a HUP record and want to inform us about that
*/
function doTellUs(acc) {
    var oDp = new RemoteDataProvider("/Traces/mailer.cgi?");
    oDp.onSuccess = x_Success;
    oDp.onError = x_Error;
    var elMsg = document.getElementById("msg-" + acc);
    var elForm = document.getElementById("form-" + acc);
    
    function x_Success(obj) {
        console.info(obj.responseText);
        if (!elMsg) return;
        if (obj.responseText == "ok") {
            elMsg.innerHTML = "Your message has been sent";
            elForm.innerHTML = "";
        }
    }
    function x_Error(obj) {
        console.info(obj.responseText);
        elMsg.innerHTML = "Error occured. Please send e-mail to custserv@nlm.nih.gov";
        elForm.innerHTML = "";
    }
    
    var oData = new Date();
    var t = oData.getTime();
    utils.createCookie("mail", t, false, '/Traces');

    var sEmail; 
    var b = false;
    sMessage = "Location: " + document.location.href + "\nAccession: " + acc;
    var a = ["pmid", "url", "citation", "email"];
    for (var i = 0; i < a.length; ++i) {
        var el = document.getElementById(a[i] + "-" + acc);
        if (!el) continue;
        el = el.value.trim();
        if (a[i] == "email") {
            if (el == "") sEmail = "nobody@ncbi.nlm.nih.gov"; // JA-3422
            else sEmail = el;
        }
        if (el) {
            sMessage += "\n" + a[i] + ": " + el;
            b = true;
        }
    }

    if (b) {
        var sRequest = "cmd=mail&subject=Hupped+" + acc + "+referenced+by+article&mail=sra_portal"
            + "&message=" + encodeURIComponent(sMessage)
            + "&email=" + sEmail
            + "&t=" + t;
        
        oDp.Get(sRequest);
    }
}
/*
 *******************************************************************************
 start global components
*/
 


utils.addEvent(window, "load", function() {
    oEls = document.querySelectorAll(".expand");
    if (oEls) for (var i = 0; i < oEls.length; i++) new Expander(oEls[i]);
    
    x = document.querySelectorAll(".taxon_id");
    for (var i = 0; i < x.length; ++i) {
        var X = new Taxonomy();
        var obj = X.Lookup(x[i]);
    }
    
});


;
jQuery(function(){
    jQuery('#rss_createfeed').bind('click',createRssFeed);
    function createRssFeed (e){
        e.preventDefault();
        var oThis = jQuery(this);
	   	var args = {
            'QueryKey': oThis.data('qk'),
            'Db': oThis.data('db'),
            'RssFeedName': jQuery('#rss_name').val(),
            'RssFeedLimit': jQuery('#rss_results').val(),
            'HID': oThis.data('hid')
        };
        Portal.$send('CreateRssFeed',args);
    }  
});

;
(function($){

    $(function() {    

        var theSearchInput = $("#term");
        var originalTerm = $.trim(theSearchInput.val());
        var theForm = jQuery("form").has(theSearchInput);
        var dbNode = theForm.find("#database");
        var currDb = dbNode.val();
        var sbConfig = {};
        try{
            sbConfig = eval("({" + theSearchInput.data("sbconfig") + "})");
        }catch(e){}
        var defaultSubmit =  sbConfig.ds == "yes";
        var searched = false;
        var dbChanged = null; //since db.change is triggered as a work around for JSL-2067 
        var searchModified = false; //this is used to allow searching when something esle changed on the page with out the term changing
    
        if(!$.ncbi)
            $.extend($,{ncbi:{}});
        if(!$.ncbi.searchbar)
            $.extend($.ncbi,{searchbar:{}});
            
        $.extend($.ncbi.searchbar,
            (function(){
                //*****************private ******************/
               function doSearchPing() {
                   try{
                    var cVals = ncbi.sg.getInstance()._cachedVals;
                    var searchDetails = {}
                    searchDetails["jsEvent"] = "search";
                    var app = cVals["ncbi_app"];
                    var db = cVals["ncbi_db"];
                    var pd = cVals["ncbi_pdid"];
                    var pc = cVals["ncbi_pcid"];
                    var sel = dbNode[0];
                    var searchDB = sel.options[sel.selectedIndex].value;
                    var searchText = theSearchInput[0].value;
                    if( app ){ searchDetails["ncbi_app"] = app.value; }
                    if( db ){ searchDetails["ncbi_db"] = db.value; }
                    if( pd ){ searchDetails["ncbi_pdid"] = pd.value; }
                    if( pc ){ searchDetails["ncbi_pcid"] = pc.value; }
                    if( searchDB ){ searchDetails["searchdb"] = searchDB;}
                    if( searchText ){ searchDetails["searchtext"] = searchText;}
                    ncbi.sg.ping( searchDetails );
                   }catch(e){
                       console.log(e);
                   }
                }
                function getSearchUrl(term){
                    var url = "";
                    if (typeof(NCBISearchBar_customSearchUrl) == "function") 
                            url = NCBISearchBar_customSearchUrl();
                    if (!url) {
                        var searchURI = dbNode.find("option:selected").data("search_uri");
                        url = searchURI ?  searchURI.replace('$',term) : 
                             "/" + dbNode.val() + "/" + ( term !="" ? "?term=" + term : "");
                        }
                    return url;
                }
            
                return {
                    //*****************exposed attributes and functions ******************/
                    'theSearchInput':theSearchInput,
                    'theForm':theForm,
                    'dbNode':dbNode,
                    'searched':searched,
                    'setSearchModified':function() { searchModified = true; },
                    'setSearchUnmodified':function() { searchModified = false; },
                    'searchModified':function(){return searchModified;},
                    'doSearch':function(e){
                           e.stopPropagation();
                           e.preventDefault();
                           //checking for the searched flag is necessary because the autocompelete control fires on enter key, the form submit also fires on enter key
                           if(searched == false){
                               searched = true;
                               theForm.find('input[type="hidden"][name^="p$"]').attr('disabled', 'disabled');
                               //$("input[name]").not(jQuery(".search_form *")).attr('disabled', 'disabled');
                               if (defaultSubmit)
                                   $.ncbi.searchbar.doSearchPing();
                               else {
                                   var term = $.trim(theSearchInput.val());
                                   if (dbChanged || searchModified || term !== originalTerm){
                                       $.ncbi.searchbar.doSearchPing();
                                       var searchUrl = $.ncbi.searchbar.getSearchUrl(encodeURIComponent(term).replace(/%20/g,'+'));
                                       var doPost = (term.length  > 2000) ? true : false; 
                                       if (doPost){
                                           if (e.data.usepjs){
                                               Portal.$send('PostFrom',{"theForm":theForm,"term":term,"targetUrl":searchUrl.replace(/\?.*/,'')});
                                           }
                                           else{
                                               theForm.attr('action',searchUrl.replace(/\?.*/,''));
                                               theForm.attr('method','post');
                                           }
                                       }
                                       else {
                                           window.location = searchUrl;
                                       }
                                   }
                                   else{ //if (term !== originalTerm){
                                       searched = false;
                                   }
                               }
                           }
                    },
                    'onDbChange':function(e){
                         if (dbChanged === null)
                             dbChanged = false;
                         else
                             dbChanged = true;
                         var optionSel = $(e.target).find("option:selected");
                         var dict = optionSel.data("ac_dict");
                         if (dict){
                             //theSearchInput.ncbiautocomplete("option","isEnabled",true).ncbiautocomplete("option","dictionary",dict);
                             theSearchInput.ncbiautocomplete().ncbiautocomplete({
                                    isEnabled: true,
                                    dictionary: dict
                                });
                             theSearchInput.attr("title","Search " + optionSel.text() + ". Use up and down arrows to choose an item from the autocomplete.");
                         }
                         else{
                           theSearchInput.ncbiautocomplete().ncbiautocomplete("turnOff",true);
                           theSearchInput.attr("title", "Search " + optionSel.text());
                         }
                         if (defaultSubmit)
                            theForm.attr('action','/' + dbNode.val() + '/');  
                    },
                    'doSearchPing':function(){
                        doSearchPing();
                    },
                    'getSearchUrl':function(term){
                        return getSearchUrl(term);
                    }
                    
                };//end of return 
             })() //end of the self executing anon
        );//end of $.extend($.ncbi.searchbar
    
         function initSearchBar(usepjs){
            //enable the controls for the back button
            theForm.find('input[type="hidden"][name^="p$"]').removeAttr('disabled');
             if (usepjs)
                 portalSearchBar();
         }
         
        
    
        function portalSearchBar(){
            
            Portal.Portlet.NcbiSearchBar = Portal.Portlet.extend ({
                init:function(path,name,notifier){
                    this.base (path, name, notifier);
                },
                send:{
                    "Cmd":null,
                    "Term":null
                },
                "listen":{
                    "PostFrom":function(sMessage,oData,sSrc){
                        this.postForm(oData.theForm,oData.term,oData.targetUrl);
                    }
                },
                "postForm":function(theForm,term,targetUrl){
                       //console.log('targetUrl = ' + targetUrl);
                       theForm.attr('action',targetUrl);
                       theForm.attr('method','post');
                       this.send.Cmd({
                            'cmd' : 'Go'
                        });
                           this.send.Term({
                            'term' : term
                        });
                        Portal.requestSubmit();
                },
                'getPortletPath':function(){
                    return this.realpath + '.Entrez_SearchBar';
                }
            });
    
        }//portalSearchBar
        


         //portal javascript is required to make a POST when the rest of the app uses portal forms 
         var usepjs = sbConfig.pjs == "yes"; 
         //console.log('sbConfig',sbConfig);
         initSearchBar(usepjs);
         
         dbNode.on("change",$.ncbi.searchbar.onDbChange);
        
        theForm.on("submit",{'usepjs':usepjs},$.ncbi.searchbar.doSearch);
        theSearchInput.on("ncbiautocompleteenter ncbiautocompleteoptionclick", function(){theForm.submit();});
        //a work around for JSL-2067
        dbNode.trigger("change");
        //iOS 8.02 changed behavior on autofocus, should probably check other mobile devices too
        if (sbConfig.afs == "yes" && !/(iPad|iPhone|iPod)/g.test(navigator.userAgent) ){ 
            window.setTimeout(function(){
                try{
                	var x = window.scrollX, y = window.scrollY; // EZ-8676
                	
                    var size= originalTerm.length;
                    if (size == 0 || /\s$/.test(originalTerm))
                        theSearchInput.focus()[0].setSelectionRange(size, size);
                    else
                        theSearchInput.focus().val(originalTerm + " ")[0].setSelectionRange(size+1, size+1);
                        
                    window.scrollTo(x, y);
                }
                catch(e){} //setSelectionRange not defined in IE8
            },1);
        }
        
        //set the query changed flag true after a few seconds, still prevents scripted clicking or stuck enter key
        window.setTimeout(function(){$.ncbi.searchbar.setSearchModified();},2000);
         
     });//End of DOM Ready

})(jQuery);

/*
a call back for the 'Turn off' link at the bottom of the auto complete list
*/
function NcbiSearchBarAutoComplCtrl(){
    jQuery("#term").ncbiautocomplete("turnOff",true);
    if (typeof(NcbiSearchBarSaveAutoCompState) == 'function')
        NcbiSearchBarSaveAutoCompState();
 }

 



;
jQuery(function () {
    Portal.Portlet.Entrez_SearchBar = Portal.Portlet.NcbiSearchBar.extend ({
        init:function(path,name,notifier){
            this.base (path, name, notifier);
            var oThis = this;
            jQuery("#database").on("change", function(){
                oThis.send.DbChanged({'db' : this.value});
            });
        },
        send:{
            "Cmd":null,
            "Term":null,
            "DbChanged":null
        },
        'listen':{
            "PostFrom":function(sMessage,oData,sSrc){
        	    this.postForm(oData.theForm,oData.term,oData.targetUrl);
        	    },
            "ChangeAutoCompleteState": function(sMessage, oData, sSrc) {
        	    this.ChangeAutoCompleteState(sMessage, oData, sSrc);
                },
            'CreateRssFeed':function(sMessage,oData,sSrc){
                this.createRssFeed(sMessage,oData,sSrc);
            },
            'AppendTerm': function(sMessage, oData, sSrc) {
    		    this.ProcessAppendTerm(sMessage, oData, sSrc);
    		},
    		// to allow any other portlet to clear term if needed  
    		'ClearSearchBarTerm': function(sMessage, oData, sSrc) {
    			jQuery("#term").val("");
    		},
    		// request current search bar term to be broadcast  
    		'SendSearchBarTerm': function(sMessage, oData, sSrc) {
    			this.send.Term({'term' : jQuery("#term").val()});
    		}
        },
        'createRssFeed':function(sMessage,oData,sSrc){
            
            var site = document.forms[0]['p$st'].value;
    	   	var portletPath = this.getPortletPath();
    	   	
            try{
                var resp = xmlHttpCall(site, portletPath, 'CreateRssFeed', oData, receiveRss, {}, this);
            }
            catch (err){
                alert ('Could not create RSS feed.');
            }
            function receiveRss(responseObject, userArgs) {
        	    try{
            	    //Handle timeouts 
            	    if(responseObject.status == 408){
            	        //display an error indicating a server timeout
            	        alert('RSS feed creation timed out.');
            	    }
            	    
            	    // deserialize the string with the JSON object 
            	    var response = '(' + responseObject.responseText + ')'; 
            	    var JSONobject = eval(response);
            	    // display link to feed
            	    jQuery('#rss_menu').html(JSONobject.Output,true);
            	    //jQuery('#rss_dropdown a.jig-ncbipopper').trigger('click');
            	    jQuery('#rss_dropdown a.jig-ncbipopper').ncbipopper('open');
            	    //document.getElementById('rss_menu').innerHTML = JSONobject.Output;
                }
                catch(e){
                    alert('RSS unavailable.');
                }
            }
                
        },
        'getPortletPath':function(){
            return this.realpath + '.Entrez_SearchBar';
        },
        "ChangeAutoCompleteState": function(sMessage, oData, sSrc){
            var site = document.forms[0]['p$st'].value;
            var resp = xmlHttpCall(site, this.getPortletPath(), "ChangeAutoCompleteState", {"ShowAutoComplete": 'false'}, function(){}, {}, this);
        },
        "ProcessAppendTerm" : function(sMessage, oData, sSrc){
            var theInput = jQuery("#term");
    	    var newTerm = theInput.val();
    	    if (newTerm != '' && oData.op != ''){
    	        newTerm = '(' + newTerm + ') ' + oData.op + ' ';
    	    }
    	    newTerm += oData.term;
    	    theInput.val(newTerm); 
    	    
    	    theInput.focus();
    	}
    }); //end of Portlet.extend
}); //end of jQuery ready

function NcbiSearchBarSaveAutoCompState(){
    Portal.$send('ChangeAutoCompleteState');
}


;
Portal.Portlet.Entrez_Facets = Portal.Portlet.extend ({
  
	init: function (path, name, notifier) 
	{ 
		this.base (path, name, notifier);
		var jFacetObj = jQuery(".facet_cont");
		if (jFacetObj[0]){
    		jFacetObj.find('.facet a').live('click',{'thisObj':this},this.filterClicked);
    		jFacetObj.find('.facet_more_apply').live('click',{'thisObj':this},this.facetMoreApplyClicked);
    		jFacetObj.find('.facet_tools a.jig-ncbipopper').bind('ncbipopperopen',{'thisObj':this},this.onMoreFilterGroups);
    		jFacetObj.find('#filter_groups_apply').bind('click',{'thisObj':this},this.filterGroupsApplyClicked);
    		jFacetObj.find('.btn_date_apply').live('click',{'thisObj':this},this.dateRangeApplyClicked);
    		jFacetObj.find('.btn_date_clear').live('click',{'thisObj':this},this.dateRangeClearClicked);
    		jFacetObj.find('.btn_range_apply').live('click',{'thisObj':this},this.rangeApplyClicked);
    		jFacetObj.find('.btn_range_clear').live('click',{'thisObj':this},this.rangeClearClicked);
    		jFacetObj.find('#facet_fields_apply').live('click',{'thisObj':this},this.facetFieldsApplyClicked);
    		
    		jFacetObj.find('.facet .more a').live('ncbipopperopen',{'thisObj':this},this.onMoreFiltersOpen);
    		jFacetObj.find('.facets_dialog').live('keypress',{'thisObj':this},this.facetDialogKeyPress);
    		jFacetObj.find('.input_date_ym').live('blur',this.autoFillDateInputs);
    		jQuery('#reset_from_message_res').live('click',{'thisObj':this},this.resetFromMessageRes);
    		
    		this.DefaultShownFacetGroups = jFacetObj.data('default_grps').split(',');
    		
    		jFacetObj.find("input[type=checkbox]").live("change",function(e){
    		   ncbi.sg.ping( this, e, "additionalFilters", { "action" : this.checked ? "checked" : "unchecked" } );
    		});
    		
    		jFacetObj.find(".of_sel_inp").live("ncbiautocompleteoptionclick", //ncbiautocompleteenter results in multiple events
    		    {'thisObj':this},this.openFieldSelected).live("keypress",{'thisObj':this},this.openFieldKeyPress);  
    		jFacetObj.find("ul.facet li.of_sel button.of_add").live("click",{'thisObj':this},this.openFieldAddClicked);
    		jFacetObj.find(".of_sel_inp").live("keyup ncbiautocompleteoptionclick input",{'thisObj':this},this.openFieldChanged);
    		
    		this.jFacetObj = jFacetObj;
    	}
		
		jQuery('#reset_from_message').on('click',{'thisObj':this},this.resetFromMessage);
		
	},
	'send':{
	    'Cmd':null,
	    'SendSearchBarTerm': null,
	    'SetTimelineFilter':null,
	    'QueryKey':null,
	    'LinkName':null,
	    'IdsFromResult':null
	},
	'listen':{
	    'FacetFilterSet':function(sMessage,oData,sSrc){
		    this.handleFacetFilterSet(oData.FacetsUrlFrag,oData.BMFacets);
		},
		'FacetFiltersCleared':function(sMessage,oData,sSrc){
		    this.handleFacetFiltersCleared();
		}
	},
	'DefaultShownFacetGroups':[],
	'jFacetObj':null,
	'filterClicked':function(e){
	    e.preventDefault();
	    var oThis = jQuery(this);
	    var facetUl = oThis.closest("ul.facet");
	    var filter_id = facetUl.data('filter_id'),value_id = oThis.data('value_id');
	    var check_on = ! oThis.parent().hasClass("selected");
	    if (value_id == 'reset'  )
	        Portal.$send('FacetFilterSet',{'FacetsUrlFrag': 'fcl=all'});
	    else if (value_id == 'fetch_more'  ){
	        if (!oThis.hasClass("jig-ncbipopper"))
	            e.data.thisObj.FetchMoreOptions(filter_id,oThis);
	    }
	    else if (value_id == 'fetch_more_exp')
	        e.data.thisObj.ShowAllFacetsToggle(e);
	    else if (filter_id == 'field_search' ){
	        if (!oThis.hasClass("jig-ncbipopper"))
	            e.data.thisObj.removeFieldSelection();
	    }
	    else if (oThis.parent().hasClass('of_sel'))
	        return;
	    else if (facetUl.data('of')=='yes' && oThis.parent().hasClass('of_fil_val')){
	        if (check_on)
	            e.data.thisObj.applyOpenField(oThis,filter_id);
	        else
	            e.data.thisObj.removeOpenField(oThis,filter_id);
	    }
	    else if (facetUl.data('of')=='yes' && !oThis.parent().hasClass('fil_val'))
	        e.data.thisObj.removeOpenField(oThis,filter_id);
	        
	    else if (facetUl.data('ss')=='yes')
	        e.data.thisObj.handleFilterSelection({'filterId':filter_id.toString(),'valueId':value_id.toString(),'checkOn':check_on,'replaceAll':true});
	    else if ((filter_id || value_id) && !oThis.hasClass("jig-ncbipopper") && !oThis.hasClass("facet_more_cancel") )
    	    e.data.thisObj.handleFilterSelection({'filterId':filter_id.toString(),'valueId':value_id.toString(),'checkOn':check_on,
    	        'dateSearch':facetUl.data('ds')=='yes','rangeSearch':facetUl.data('rs')=='yes'});
    	
        
	},
    'handleFilterSelection':function(opts){
	    var defOpts = {'filterId':undefined,'valueId':undefined,'checkOn':undefined,'replaceAll':undefined,'dateSearch':undefined,'rangeSearch':undefined};
	    opts = jQuery.extend(defOpts,opts);
	    
	    //when replaceAll is true, all values in that filter group are replaced, used for single select groups
	    //valueId == ''  means clear that group 
	    //var currFilterString = window.location.search.match(/filters=([^&]*)/);
	    var currFilterString = this.getValue('FacetsUrlFrag').match(/filters=([^&]*)/);
	    //var currFilterVals = currFilterString && currFilterString[1] ? currFilterString[1].split(';') : [];
	    var currFilterVals = currFilterString ? currFilterString[1].split(';') : [];
	    var possibleVals = [];
	    var facetGrpUl = this.jFacetObj.find('ul[data-filter_id = "' + opts.filterId + '"]');
	    facetGrpUl.find('li.fil_val a').each(function(){
	        var possIdVal = jQuery(this).data('value_id');
	        if (possIdVal)
	            possibleVals.push(possIdVal.toString());
	        });
	    currFilterVals = this.customFilterRead(currFilterVals,possibleVals,opts.filterId,opts.dateSearch,opts.rangeSearch);
	    
	    function removeValues(valuesArr) {
	        jQuery(valuesArr).each(function(ind,val){
	            var indexInCurr = jQuery.inArray(val,currFilterVals);
	            if (indexInCurr != -1)
	                 currFilterVals.splice(indexInCurr,1);
	        });
	    }
	    function addValues(valuesArr) {
	        jQuery(valuesArr).each(function(ind,val){
	             var indexInCurr = jQuery.inArray(val,currFilterVals);
	             if (indexInCurr == -1)
	                 currFilterVals.push(val);
	        });
	    }
	    
	    if (opts.replaceAll == true && opts.checkOn){ //single select
	        removeValues(possibleVals);
	        addValues(opts.valueId.split(';'));
	    }
	    else if (opts.valueId == ''){
	        removeValues(possibleVals);
	    }
	    else if (opts.checkOn){
	        addValues(opts.valueId.split(';'));
	    }
	    else if (!opts.checkOn){
	        removeValues(opts.valueId.split(';'));
	    }
	    var bmFacets = '';
	    if (facetGrpUl.data('bm') == 'yes' && !(opts.checkOn != true && facetGrpUl.find('li.selected').size() == 1) ){
	        bmFacets = 'bmf=' + facetGrpUl.data('filter_id') + ':' +
	            jQuery.makeArray(facetGrpUl.find('li.fil_val a').map(function(){return (jQuery(this).data('value_id'))})).join(';');
	    }
	    
	    Portal.$send('FacetFilterSet',{'FacetsUrlFrag':this.getNewUrlFrag(currFilterVals.join(';')),'BMFacets':bmFacets});
        
	},	
	'customFilterRead':function(currFilterVals,possibleVals,filterId,datesearch,rangesearch){
	    //if there is db specific filter reading override this
	    if(datesearch == true){ 
	        var rg = new RegExp(filterId + '_' + '\\d{4}\/\\d{2}\/\\d{2}_\\d{4}\/\\d{2}\/\\d{2}');
	        //for (var ind in currFilterVals){
	        for(var ind=0; ind<currFilterVals.length; ind++){
	            if (rg.exec(currFilterVals[ind]) ||
	                jQuery.inArray(currFilterVals[ind],possibleVals) != -1 ){
	                currFilterVals.splice(ind,1);
	            }
	        }
	    }
	    else if (rangesearch == true){
	        var rg = new RegExp(filterId + '_[^_]+_[^_]+');
	        for(var ind=0; ind<currFilterVals.length; ind++){
	            if (rg.exec(currFilterVals[ind]) ||
	                jQuery.inArray(currFilterVals[ind],possibleVals) != -1 ){
	                currFilterVals.splice(ind,1);
	            }
	        }
	    }
	    return currFilterVals;
	},
	'getNewUrl':function(filters,fcl,allowEmptyTerm){
	    var currUrl = window.location.pathname + window.location.search ;
        currUrl = this.replaceUrlParam(currUrl, 'filters', filters);  
        currUrl = this.replaceUrlParam(currUrl,'fcl', fcl); 
        currUrl = this.replaceUrlParam(currUrl,'querykey','');
        currUrl = this.replaceUrlParam(currUrl,'cmd','');
        currUrl = this.addTermToUrl(currUrl,allowEmptyTerm);
        //currUrl = this.appendUrlHash(currUrl);
        return currUrl;
	},
	'addTermToUrl':function(currUrl,allowEmptyTerm){
/*	    if (!currUrl.match(/term=*\/)){
	        //currUrl = this.replaceUrlParam(currUrl,'term',this.jFacetObj.data('term'));
	    } */
	    var term = jQuery.trim(jQuery("#search_term").val());
	    if (allowEmptyTerm != true)
	        term = term == '' ? 'all[sb]' : term;
	    currUrl = this.replaceUrlParam(currUrl,'term',term);
	    return currUrl;
	},
	'replaceUrlParam':function(currUrl,paramName,paramVal,allowEmpty){
	    paramVal = paramVal ? paramVal : '';
        if (paramVal != '' || allowEmpty)
            if (currUrl.indexOf(paramName + '=') == -1)
                currUrl = currUrl + (currUrl.indexOf('?') != -1 ? '&' : '?') + paramName + '=' + paramVal;
            else
                currUrl = currUrl.replace(new RegExp(paramName + '=[^&]*'), paramName + '=' + paramVal);
         else
             if (currUrl.match(new RegExp('&' + paramName + '=[^&]*')))
                 currUrl = currUrl.replace(new RegExp('&' + paramName + '=[^&]*'),'');
             else if (currUrl.match(new RegExp(paramName + '=[^&]*&')))
                 currUrl = currUrl.replace(new RegExp(paramName + '=[^&]*&'),'');
             else
                 currUrl = currUrl.replace(new RegExp(paramName + '=[^&]*'),'');
         return currUrl;
	},
	'getNewUrlFrag':function(filters,fcl){
	    var currUrl = this.getValue('FacetsUrlFrag');
        currUrl = this.replaceParamFrag(currUrl, 'filters', filters);
        currUrl = this.replaceUrlParam(currUrl,'fcl', fcl); 
        return currUrl;
	},
	'replaceParamFrag':function(currUrl,paramName,paramVal){//TO-DO ... poorly named, refactor
          //currUrl = currUrl.replace(new RegExp(paramName + '=[^;]*'), paramName + '=' + paramVal);
          currUrl = 'filters=' + paramVal;
          return currUrl;
	},
	'replaceUrlParamFrag':function(origFrag,paramName,paramVal,delim){ 
	    delim = delim || ';';
	    if (paramVal != '')
            if (origFrag.indexOf(paramName + '=') == -1)
                return  origFrag == '' ? paramName + '=' + paramVal : origFrag + delim + paramName + '=' + paramVal ;
            else
                return origFrag.replace(new RegExp(paramName + '=.[^' + delim + ']*'), paramName + '=' + paramVal);
         else
             if (origFrag.match(new RegExp(delim + paramName + '=.[^' + delim + ']*')))
                 return origFrag.replace(new RegExp(delim + paramName + '=.[^' + delim + ']*'),'');
             else if (origFrag.match(new RegExp(paramName + '=.[^' + delim + ']*' + delim)))
                 return origFrag.replace(new RegExp(paramName + '=.[^' + delim + ']*' + delim),'');
             else 
                 return origFrag.replace(new RegExp(paramName + '=.[^' + delim + ']*'),'');
        
	},
	'appendUrlHash':function(urlStr){
	    var hash = window.location.hash;
        if (hash != '')
            urlStr = urlStr + "#" + hash;
        return urlStr;
	},
	'FetchMoreOptions':function(filter_id,moreNode){
	    //if the moreNode param is not null, coming from a 'more' under a category, otherwise it is adding a whole group from 'choose filters'
	    var args = {"MoreFacetsGroupId":filter_id,"MoreFacetsNewGroup":(moreNode?"":"true"),"Db":this.jFacetObj.data('db'),"Term":jQuery("#term").val()};
        var site = document.forms[0]['p$st'].value;
        // ajax call
        xmlHttpCall(site, this.getPortletPath(), "GetMoreFilters", args, this.receiveMoreFilters, {"moreNode":moreNode}, this);
	},
	'receiveMoreFilters':function(responseObject, userArgs){
        try {
            // Handle timeouts
            if (responseObject.status == 408) {
                //this.showMessage("Server currently unavailable. Please check connection and try again.","error");
                 console.warn("Server currently unavailable. Please check connection and try again.");
                return;
            }
            var resp = '(' + responseObject.responseText + ')';
            var JSONobj = eval(resp);
            var allFilters = JSONobj.all_filters;
            if (userArgs.moreNode)
                this.addMoreFiltersDialog(allFilters,userArgs.moreNode);
            else
                this.addMoreFilterGroup(allFilters);
            //TO-DO: be more specific about this scan
            jQuery.ui.jig.scan();
            
        } catch (e) {
            //this.showMessage("Server error: " + e, "error");
            console.warn("Server error: " + e);
        }
	},
	'addMoreFiltersDialog':function(allFilters,targetNode){
	    targetNode.addClass("jig-ncbipopper");
	    var popper = jQuery(targetNode.attr('href'));
	    var filterId = targetNode.closest("ul.facet").data('filter_id');
	    var selFilters = this.jFacetObj.find('ul[data-filter_id = "' + filterId + '"] li a');
	    allFilters = jQuery(allFilters);
	    selFilters.each(function(){
	        allFilters.find('li input[id = "' + jQuery(this).data('value_id') + '"]').attr('checked','checked');
	        });   
	    popper.append(allFilters);
	    jQuery.ui.jig.scan(targetNode,['ncbipopper']);
	    targetNode.ncbipopper('open');
	},
	'getPortletPath': function(){
        return this.realname;
    },
    'facetMoreApplyClicked':function(e){
        e.preventDefault();
        var self = jQuery(e.target);
        if (self.find('span').text() == 'Add'){
            e.data.thisObj.addOpenFieldValue(self.closest('ul.facet'));
            return;            
        }
        var facetGroup = self.closest('ul.facet');
        var groupId = facetGroup.data('filter_id');
        var selFilters = jQuery('#' + groupId + '_more').find('li input').filter('input:checked');
        var filtersInFacet = facetGroup.find('li.fil_val a');
        var ofFiltersInFacet = facetGroup.find('li.of_fil_val a');
        var addedFacets = [], removedFacets = [], newFacets = [];
        var isOpenField = facetGroup.find('.filter_grp').is('.of_grp');
        //alert(isOpenField);
        selFilters.each(function () {
            var oThis = jQuery(this);
            var filterId = oThis.data('value_id');
            var filterName = oThis.next().text();
            addedFacets.push(filterId);
            var parentValueId = oThis.parent().data('value_id');
            if( oThis.parent().data('value_id') == "of_val" && ofFiltersInFacet.filter(function(ind,el){return el.text == filterName;} ).size() == 0){
                jQuery('<li class="of_fil_val"><a data-qval="' + filterName + '" data-value_id="' + filterName + '" href="#">' + filterName + '</a></li>').insertBefore(facetGroup.find("li.more"));
            }
            else if (oThis.parent().data('value_id') != "of_val" && filtersInFacet.filter('a[data-value_id = "' + filterId + '"]').size() === 0){
                newFacets.push(filterId);
                //find the place to insert
                var insertBeforeNode;
                facetGroup.find('li.fil_val').each(function(){
                    if (jQuery(this).find('a').text() > filterName){
                        insertBeforeNode = jQuery(this);
                        return false;
                    }
                });
                if (!insertBeforeNode)
                    insertBeforeNode = facetGroup.find("li.more")
                    
                jQuery('<li class="fil_val"><a data-value_id="' + filterId + '" href="#">' + filterName + '</a></li>').insertBefore(insertBeforeNode);
            }
        });
        filtersInFacet.add(ofFiltersInFacet).each(function(){
            var oThis = jQuery(this);
            var filterId = oThis.data('value_id');
            if (selFilters.filter('input[data-value_id="' + filterId + '"]').size() === 0){
                removedFacets.push(filterId);
                facetGroup.find('li.fil_val').add(facetGroup.find('li.of_fil_val')).has('a[data-value_id="' + filterId + '"]').remove();
            }
        });
        
        ncbi.sg.ping( e.target, e, "additionalFiltersApply", {"allChecked" : addedFacets, "newChecked" : newFacets , "newUnchecked": removedFacets} );
        
        facetGroup.find('li a[data-value_id="fetch_more"]').ncbipopper('close');
        
        function arrayToXml(arr){
            var xmlStr = '<Facets><FacetGroup '  + ' id = "' + groupId + '" >';
            for(var ind=arr.length -1; ind >=0 ; ind--)
                xmlStr = xmlStr + '<Facet>' + arr[ind] + '</Facet>';
            xmlStr = xmlStr + '</FacetGroup></Facets>';
            return xmlStr;
        }
        var args = {"UserSelectedFacetsNew":arrayToXml(addedFacets),"UserDeSelectedFacetsNew":arrayToXml(removedFacets)};
        
        
        var site = document.forms[0]['p$st'].value;
        // ajax call
        xmlHttpCall(site, e.data.thisObj.getPortletPath(), "UpdateUserAddedFacets", args, function(){}, null, this);       
    },
    'onMoreFilterGroups':function(e){
        jQuery('#filter_groups_apply').data('attachedTo',e.target.id);
        
        var loadedFgIds = [],activeFgIds = [];
        e.data.thisObj.jFacetObj.find('.facet .filter_grp a.clear').each(function(){
            var filterGrp = jQuery(this).closest('ul.facet');
            var filterId = 'fg_' + filterGrp.data('filter_id');
            loadedFgIds.push(filterId);
            if (filterGrp.find('li.selected')[0])
                activeFgIds.push(filterId);
        });
        var fgChecks = jQuery('#more_filter_groups input');
        fgChecks.each(function(){
            var oThis = jQuery(this);
            var currId = oThis.attr('id');
            oThis.attr('checked',jQuery.inArray(currId,loadedFgIds) != -1);
            oThis.attr('disabled',oThis.data('always_show') == 'yes' || jQuery.inArray(currId,activeFgIds) != -1)
        });
    },
    'filterGroupsApplyClicked':function(e){
        e.preventDefault();
        var loadedFgIds = [], fgIdsAdd = [],fgIdsRemove = [],selFgIds = [],fgUserSelIds=[];
        var defaultShownFacetGroups = e.data.thisObj.DefaultShownFacetGroups;
        e.data.thisObj.jFacetObj.find('.facet .filter_grp a.clear').each(function(){
            loadedFgIds.push('fg_' + jQuery(this).closest('ul.facet').data('filter_id'));
        });
        e.data.thisObj.jFacetObj.find('#more_filter_groups input').filter('input:checked').each(function(){
            selFgIds.push(jQuery(this).attr('id'));
        });
        var last = selFgIds.length;
        for (var ind =0; ind <last; ind++  ){
            if(jQuery.inArray(selFgIds[ind],loadedFgIds) == -1)
                fgIdsAdd.push(selFgIds[ind].substring(3));
            if(jQuery.inArray(selFgIds[ind],defaultShownFacetGroups) == -1)
                fgUserSelIds.push(selFgIds[ind].substring(3));
        }
        last = loadedFgIds.length;
        for (var ind =0; ind <last; ind++  )
            if (jQuery.inArray(loadedFgIds[ind],selFgIds) == -1)
                fgIdsRemove.push(loadedFgIds[ind].substring(3));
        
        e.data.thisObj.updateFiltersShown(fgIdsAdd,fgIdsRemove,fgUserSelIds);
        jQuery('#' + jQuery(this).data('attachedTo')).ncbipopper('close');
    },
    'updateFiltersShown':function(fgIdsAdd,fgIdsRemove,fgUserSelIds){
        var last = fgIdsRemove.length;
        for (var ind =0; ind <last; ind++  )
            this.jFacetObj.find('ul.facet[data-filter_id = ' + fgIdsRemove[ind] + ']').remove();
        last = fgIdsAdd.length -1;
        for (var ind = last; ind >= 0; ind--  )
            this.FetchMoreOptions(fgIdsAdd[ind],null);
        //update the selection on the session variables
        this.updateUserSelectionAttrs(fgUserSelIds,fgIdsRemove);
    },
    'updateUserSelectionAttrs':function(fgUserSelIds,fgIdsRemove){
        
        function arrayToXml(arr,rootTag,tag){
            var xmlStr = '<' + rootTag + '>';
            var last = arr.length;
            for(var i=0; i<last; i++)
                xmlStr = xmlStr + '<' + tag + '>' + arr[i] + '</' + tag + '>';
            xmlStr = xmlStr + '</' + rootTag + '>';
            return xmlStr;
        }
        var rootTag = 'FacetGroups',tag='FacetGroup';
        var args = {"UserSelectedFacetGroups":arrayToXml(fgUserSelIds,rootTag,tag),"UserDeSelectedFacetGroups":arrayToXml(fgIdsRemove,rootTag,tag)};
        var site = document.forms[0]['p$st'].value;
        // ajax call
        xmlHttpCall(site, this.getPortletPath(), "UpdateUserSelectedFacetGroups", args, function(){} , {}, this);
        
    },
    'addMoreFilterGroup':function(allFilters){
	    allFilters = jQuery(allFilters);
	    
	    //console.log('addMoreFilterGroup');
	    
/*	    if(!allFilters.find("ul>li")[0]){
	        alert("That wouldn't return any results");
	        return;
	    }*/
	    
	    //find the position and insert
	    var nFilterId = allFilters.data("filter_id");
	    //console.log('curr filter id ', nFilterId);
	    var nFilerLi = jQuery('#more_filter_groups input').filter(function(i,j){return jQuery(j).attr("id") == "fg_" + nFilterId;}).parent();
	    //console.log('curr li in more dialog',nFilerLi);
	    var selFacet = nFilerLi.nextAll("li").filter(function(i,j){return jQuery(j).find("input").is(':checked')})[0];
	    //var selFacet = nFilerLi.nextAll("li").filter(function(i,j){console.log('find next sel',jQuery(j),jQuery(j).find("input").is(':checked'),jQuery(j).find("input[checked]"),jQuery(j).find("input[checked]")[0]); return jQuery(j).find("input").is(':checked')})[0];
	    //console.log('sel facet after',selFacet);
	    var facetUl;
	    if (selFacet){
	        selFacet = jQuery(selFacet);
	        var facetId = selFacet.find("input").attr("id").substring(3);
	        facetUl = jQuery("ul.facet").filter(function(i,j){return jQuery(j).data("filter_id") == facetId})
	        console.log('sel facet after ul',facetUl);
	        
	    }
	    if (facetUl && facetUl[0])
	        facetUl.before(allFilters);
	    else{
	        var resetLink = jQuery('ul.facet_reset').has('li a[data-value_id="reset"]');
	        resetLink.before(allFilters);
	    }
	    
	    var moreLink = allFilters.find("li.more");
	    if (moreLink[0]){
	        moreLink.find("a").addClass("jig-ncbipopper");
	        jQuery.ui.jig.scan(moreLink,['ncbipopper'])
	    }
	    if (allFilters.find("#facet_fileds_popup")[0])
	        jQuery.ui.jig.scan(allFilters,['ncbipopper']);
	    
	    

    },
    'rangeApplyClicked':function(e){
        e.preventDefault();
        var elem = jQuery(e.target);
        var outerDiv = elem.closest('[id^=facet_range_div]');
        var valSt = outerDiv.find('[id^=facet_range_st]').val();
        var valEnd = outerDiv.find('[id^=facet_range_end]').val();
        var filterId = outerDiv.closest('ul.facet').data('filter_id');
        
        function validate(){
            var valid = true;
            try{
                var validationRE = outerDiv.data('vre') || '[^\s]+';
                var rg = new RegExp(validationRE);
                valid = valid && Boolean(rg.exec(valSt)) && Boolean(rg.exec(valEnd));
                
                //now check for value function
                var valueFun = outerDiv.data('vf');
                if (valueFun && valid){
                    valueFunEval = eval('(' + valueFun + ')');
                    if(typeof valueFunEval == 'function')
                        valid =  valueFunEval(valEnd) > valueFunEval(valSt); 
                    else{
                        var stValue = valueFun.replace('$',valSt);
                        stValue=eval('(' + stValue + ')');
                        var endValue = valueFun.replace('$',valEnd);
                        endValue = eval('(' + endValue + ')');
                        valid = endValue >= stValue;
                    }
                }
            }
            catch(e){
                alert('Check your validation regular expressions and functions in the source xml. Your user should never see this!');
                console.error(e);
                return false;
            }
            
            return valid;
        }
        
        var tryAgain = !(e.data.thisObj.validateRange(outerDiv) && validate()); 
        if (tryAgain){
	        alert('please enter a valid range');
	        return;
	    }
	    rangeValue = filterId + '_' + valSt + '_' + valEnd;
	    e.data.thisObj.handleFilterSelection({'filterId':filterId,'valueId':rangeValue,'checkOn':true,'rangeSearch':true}); 
	    outerDiv.data('attached-to').ncbipopper('close');
    },
    //this function is a callback. If you want to have extra validation of range values - override
    'validateRange':function(outerDiv){
        return true;
    },
    'dateRangeApplyClicked':function(e){
        e.preventDefault();
        var dateRange = '',dateRangeVals = [],tryAgain = false;
        
        //if (fieldSize == 4){
        var fieldSize = 4;
        //var year1 = jQuery('#facet_date_st_year');
        var outerDiv = jQuery(e.target).closest("[id^=facet_date_range_div]");
        var year1 = outerDiv.find('[id^=facet_date_st_year]');
        //var year2 = jQuery('#facet_date_end_year');
        var year2 = outerDiv.find('[id^=facet_date_end_year]');
        var year1Val = year1.ncbiplaceholder().ncbiplaceholder('value');
        var year2Val = year2.ncbiplaceholder().ncbiplaceholder('value');
        var year1Okay = year1Val.match(new RegExp('^\\d{' + fieldSize + '}$'));
        var year2Okay = year2Val.match(new RegExp('^\\d{' + fieldSize + '}$'));
        var oneYearThere = false;
        if (year1Val == '' && year2Okay){
            year1.val('0001');
            oneYearThere = true;
        }
        else if (year2Val == '' && year1Okay){
            year2.val('3000');
            oneYearThere = true;
        }
        if ( !oneYearThere  &&  !(year1Okay && year2Okay) )
            tryAgain = true;

        if (!tryAgain){
           //jQuery('#facet_date_range_div input').each(function(){
           outerDiv.find('input').each(function(){
                var oThis = jQuery(this);
                var val = oThis.ncbiplaceholder().ncbiplaceholder('value'); //.val();
                var fieldSize = oThis.attr('size');
                if(this.id.match('month')){
                    if (!val.match(new RegExp('^\\d{0,' + fieldSize + '}$')) )
                        tryAgain = true;
                    else if (val == '' )
                        val = this.id.match("end") ? '12' : '01' ;
                    else if (val.length == 1) 
                        val = '0' + val;
                    else if (Number(val) > 12)
                        tryAgain = true;
                }
                else if(this.id.match('day')){
                    if (!val.match(new RegExp('^\\d{0,' + fieldSize + '}$')) )
                        tryAgain = true;
                    else if (val == '' )
                        val = this.id.match("end") ? '31' : '01' ;
                    else if (val.length == 1) 
                        val = '0' + val;
                    else if (Number(val) > 31)
                        tryAgain = true;
                }
                dateRangeVals.push(val);
            });
        }
	    if (tryAgain){
	        alert('please enter a valid date range');
	        return;
	    }
	    var filterId = outerDiv.closest('ul.facet').data('filter_id');
	    dateRange = filterId + '_' + dateRangeVals[0] + '/' + dateRangeVals[1] + '/' + dateRangeVals[2] + '_' + dateRangeVals[3] + '/' + dateRangeVals[4] + '/' + dateRangeVals[5];
	    e.data.thisObj.handleFilterSelection({'filterId':filterId,'valueId':dateRange,'checkOn':true,'dateSearch':true});
	    outerDiv.data('attached-to').ncbipopper('close');
	},
	'facetFieldsApplyClicked':function(e){
	    e.preventDefault();
	    var val = jQuery('#facet_fileds_select').val();
	    //var currFilterString = window.location.search.match(/filters=([^&]*)/);
	    var currFilterString = e.data.thisObj.getCurrentFilterString();
	    if (currFilterString.match(/fld_.+/)){
	        currFilterString = currFilterString.replace(/fld_.[^;]+/,val);       
	    }
	    else
	        currFilterString = (currFilterString != '') ? currFilterString + ';' + val : val; 
	    Portal.$send('FacetFilterSet',{'FacetsUrlFrag':e.data.thisObj.getNewUrlFrag(currFilterString)});
	},
	'removeFieldSelection':function(){
	    //var currUrl = window.location.pathname + window.location.search ;
	    var currUrl = this.getValue('FacetsUrlFrag');
         if (currUrl.match(/;fld_.[^;]+/))
             currUrl = currUrl.replace(/;fld_.[^;]+/,'');
         else if (currUrl.match(/fld_.[^;]+;/))
             currUrl = currUrl.replace(/fld_.[^;]+;/,'');
         else if (currUrl.match(/fld_.[^;]+/))
             currUrl = currUrl.replace(/fld_.[^;]+/,''); 
         currUrl = this.getNewUrlFrag(currUrl);
         Portal.$send('FacetFilterSet',{'FacetsUrlFrag':currUrl});
         //window.location = currUrl;
	},
	'onMoreFiltersOpen':function(e){
	    var targetNode = jQuery(this);
	    var popper = jQuery(targetNode.attr('href'));
	    var filterId = targetNode.closest("ul.facet").data('filter_id');
	    var facetUl = e.data.thisObj.jFacetObj.find('ul[data-filter_id = "' + filterId + '"]');
	    var selFilters = facetUl.find('li.fil_val a');
	    selFilters = selFilters.add(facetUl.find('li.of_fil_val a'));
	    selFilters.each(function(){
	        var self = jQuery(this);
	        popper.find('li input[data-value_id = "' + jQuery(this).data('value_id') + '"]').attr('checked','checked');
	        }); 
	    var activeFilters = selFilters.filter(function(){return jQuery(this).parent().hasClass("selected");});
	    activeFilters.each(function(){
	        popper.find('li input[data-value_id = "' + jQuery(this).data('value_id') + '"]').attr('disabled','true');
	    });
	},
	'facetDialogKeyPress':function(e){
	    e = e || utils.fixEvent (window.event);
	    if ((e.keyCode || e.which) == 13){
	        e.preventDefault();
	        jQuery(this).find('button.primary-action').trigger('click');
	    }
	},
	'autoFillDateInputs':function(e){
	    var oThis = jQuery(this);
	    var outerDiv = oThis.closest('[id^=facet_date_range_div]');
	    function updateVal(jSel,value){
	        jSel.each(function(){ var oThis = jQuery(this); if (oThis.val() == '') oThis.val(value);});
	    }
	    if (oThis.val().match(new RegExp('^\\d{' + oThis.attr('size') +'}$'))){
	        var currId = oThis.attr('id');
	        if( currId.match(/^facet_date_st_year/))
	            updateVal(outerDiv.find('[id^=facet_date_st_month], [id^=facet_date_st_day]'),'01');
	        else if (currId.match(/^facet_date_st_month/))
	            updateVal(outerDiv.find('[id^=facet_date_st_day]'),'01');    
	        else if (currId.match(/^facet_date_end_year/)){
	            updateVal(outerDiv.find('[id^=facet_date_end_month]'),'12');
	            updateVal(outerDiv.find('[id^=facet_date_end_day]'),'31');
	        }
	        else if (currId.match(/^facet_date_end_month/))
	            updateVal(outerDiv.find('[id^=facet_date_end_day]'),'31'); 
	    }
	},
	'dateRangeClearClicked':function(e){
	    e.preventDefault();
	    var self = jQuery(e.target);
	    if (self.closest('ul').has('li.daterange').find('li.selected')[0])
	        e.data.thisObj.handleFilterSelection({'filterId':self.closest('ul.facet').data('filter_id'),'valueId':'','checkOn':true,'dateSearch':true});
	    else
	        self.closest('.facets_dialog').find('input').val('');
	},
	'rangeClearClicked':function(e){
	    e.preventDefault();
	    e.data.thisObj.handleFilterSelection({'filterId':jQuery(e.target).closest('ul.facet').data('filter_id'),'valueId':'','checkOn':true,'rangeSearch':true});
	},
	'resetFromMessage':function(e){
	    e.preventDefault();
	    Portal.$send('FacetFiltersCleared',{});
	},
	'resetFromMessageRes':function(e){
	    e.preventDefault();
	    Portal.$send('FacetFilterSet',{'FacetsUrlFrag': 'fcl=all'});
	},
	'getFacetSearchData':function(){
	    var sd = {};
	    try{
	        sd = eval('({' + this.jFacetObj.data('sd') + '})');
	    }catch(e){}
	    return sd;
	},
	'handleFacetFilterSet':function(facetsUrlFrag,bMFacets){
	    var sd = this.getFacetSearchData();
	    this.setValue('FacetsUrlFrag',facetsUrlFrag);
	    this.setValue('FacetSubmitted','true');
	    this.setValue('BMFacets',bMFacets);
	    this.send.SetTimelineFilter({'TimelineYear':''});
	    if(sd.extra){
	        this.handleExtraSD(sd.extra);
	    }
	    else if (sd.op == 'search'){
	        this.send.SendSearchBarTerm();
	        this.send.Cmd({'cmd':'search'});    
	    }
	    else if (sd.op == 'link' && sd.linkname && (sd.qk || sd.idsfromresult) ){
	        this.send.LinkName({'linkname':sd.linkname});
	        this.send.QueryKey({'qk':sd.qk});
	        this.send.IdsFromResult({'IdsFromResult':sd.idsfromresult});
	        this.send.Cmd({'cmd':'Link'});    
	    }
	    else{
	        this.send.Cmd({'cmd':'HistorySearch'});
	        this.send.QueryKey({'qk':sd.qk});
	    }

	    Portal.requestSubmit();
	},
	'handleExtraSD':function(extraSD){
	    alert('Please implement the function handleExtraSD');    
	},
	'handleFacetFiltersCleared':function(){
	    this.send.Cmd({'cmd': 'removefacets'});
		Portal.requestSubmit();
	},
	'openFieldSelected':function(e){
	    e.preventDefault();
        e.data.thisObj.addOpenFieldValue(jQuery(e.target).closest('ul.facet'));
	},
	'openFieldAddClicked':function(e){
	    e.preventDefault();
	    e.data.thisObj.addOpenFieldValue(jQuery(e.target).closest('ul.facet'));
	},
	'openFieldKeyPress':function(e){
	    //e.data.thisObj.openFieldChanged(e);
	    e = e || utils.fixEvent (window.event);
	    if ((e.keyCode || e.which) == 13){
	        e.preventDefault();
	        e.data.thisObj.addOpenFieldValue(jQuery(e.target).closest('ul.facet'));
	    }
	},
	'openFieldChanged':function(e){
	    var self = jQuery(this);
	    var applyBtn = self.closest('.facets_dialog').find('.facet_more_apply');
	    if(self.val() == ''){
	        applyBtn.find('span').text('Show');
	    }
	    else{
	        applyBtn.find('span').text('Add');
	    }
	},
	'checkSelOnlyOpenField':function(input,showAlert){
	      showAlert = showAlert || 'yes';
	      var isInDict = false;
	      var inputText = input.val().toLowerCase();
	      if(input.data('so') == 'yes'){
	          var jigOpts = input.data('jigconfig').match(/dictionary:'(\w+)'.*/);
	          var dict = jigOpts ? jigOpts[1] : null;
	          jigOpts = input.data('jigconfig').match(/localData:(')?([^,]*)(')?/);
	          var localDict = jigOpts ? jigOpts[2] : null;
	          if (dict){
	              var ajaxCall = jQuery.ajax({
	                  url:'/portal/utils/autocomp.fcgi?dict=' + dict + '&q=' + inputText,
	                  async:false,
	                  dataType:'json'
	              }).always(function(data){
	                  isInDict = eval(data.responseText);
	                  //the handling function with local scope only
	                  function NSuggest_CreateData(q,matches,count){
	                      var rg = new RegExp('^' + inputText + '(@.*)?$','i');
	                      return jQuery.grep(matches,function(e,i){
	                          return rg.exec(e);
	                          }).length > 0;
	                  }
	              });
        	      if (!isInDict && showAlert == 'yes')
	                  alert('Please select one of the valid values');
	              return isInDict;
	           }
	           else if (localDict){
	               var localDictSplitted = localDict.split('.');
	               var localDictVar = null;
	               for(var i=0; i<localDictSplitted.length; i++){
	                   if (localDictVar == null)
	                       localDictVar = window[localDictSplitted[i]];
	                    else
	                        localDictVar = localDictVar[localDictSplitted[i]];
	               }
	               var rg = new RegExp('^' + inputText + '$', 'i');
	               jQuery.each(localDictVar,function(ind,val){
	                   if (val.match(rg))
	                       isInDict = true;
	               });
                 if (!isInDict && showAlert == 'yes')
	                  alert('Please select one of the valid values');
	              return isInDict;
	           }
	           else
	               return true;
	       }
	       else
	           return true;
	},
	'addOpenFieldValue':function(facetUl){
	    var inputBox = facetUl.find(".of_sel_inp");
	    var newVal = inputBox.val();
	    if(newVal){
            if(!this.checkSelOnlyOpenField(inputBox)){
	            return;
	        }
	        var listUl = facetUl.find('.facets_dialog ul.facet_more');
	        if (listUl.find('li').has('input[data-value_id="' + newVal +'"]').size() == 0 ){ 
                inputBox.val('');
    	        var elId = 'ofv_' + newVal;
    	        listUl.append('<li data-value_id="of_val"><input type="checkbox" id="'+ elId +'" checked="checked" data-value_id="' + newVal + '" ><label for="'+elId+'">' + newVal +'</label></li>');
    	        inputBox.focus();
	        }
	        else{
	            alert('Already added');
	            inputBox.focus();
	        }
	    }
	    else{
	        facetUl.find('.facet_more_apply').trigger('click');
	    }
	},
	'getCurrentFilterString':function(){
	    var currFilterString = this.getValue('FacetsUrlFrag').match(/filters=([^&]*)/);
	    return currFilterString ? currFilterString[1] : ''; 
	},
	'applyOpenField':function(elem,filterId){
        var currFilterString = this.getCurrentFilterString();
        var paramVal = '';
        var newVal = elem.data('value_id');
        var dupl = false;
        var facetUl = elem.closest('ul.facet');
        facetUl.find('li.selected').not(".fil_val").each(function(){
            var currVal = jQuery(this).find('a').data('qval');
            if (newVal.match(new RegExp('^' + currVal + '$','i')))
                dupl = true;
            paramVal = paramVal + ( paramVal == '' ? '' : ':' ) + currVal ;
        });
        if (dupl)
            return;
        paramVal = paramVal == '' ? newVal : paramVal + ':' + newVal;
        currFilterString = this.replaceUrlParamFrag(currFilterString,'of_' + filterId,paramVal,';');
    	    
        
        var bmFacets = '';
        var facetUl = elem.closest('ul.facet');
        if (facetUl.data('bm') == 'yes'){
            bmFacets = 'bmf=' + facetUl.data('filter_id') + ':' +
                jQuery.makeArray(facetUl.find('li a').map(function(){return (jQuery(this).data('value_id'))})).join(';');
        }
        	    
        Portal.$send('FacetFilterSet',{'FacetsUrlFrag':this.getNewUrlFrag(currFilterString),'BMFacets':bmFacets});
	},
	'removeOpenField':function(elem,filterId){
	    var currFilterString = this.getCurrentFilterString();
	    var valueId = elem.data('value_id');

            
        var toReplace = currFilterString.match(new RegExp('of_' + filterId + '=(.[^;]*)'));
        toReplace = toReplace ? toReplace[1] : '';
        var replaceWith = '';
        if (valueId != ''){
            var toRemove = elem.data('qval');
            replaceWith = toReplace;
            var rg;
            rg = new RegExp(':' + toRemove);
            if(rg.exec(replaceWith))
                replaceWith = replaceWith.replace(rg,'');
            else{
                rg = new RegExp(toRemove + ':');
                if (rg.exec(replaceWith))
                    replaceWith = replaceWith.replace(rg,'');
                else{
                    replaceWith = replaceWith.replace(new RegExp(toRemove),'');
                }
            }
            
            
        }
        currFilterString = this.replaceUrlParamFrag(currFilterString,'of_' + filterId,replaceWith,';')
        this.setValue('FacetsUrlFrag',"filters=" + currFilterString);
        this.handleFilterSelection({'filterId':filterId,'valueId':valueId,'checkOn':true});
	},
	'ShowAllFacetsToggle':function(e){
	    var elem = jQuery(e.target);
	    if (elem.hasClass('fetch_more_exp')){
	        elem.removeClass('fetch_more_exp');
	        elem.addClass('fetch_more_exp_less');
	        if (isNaN(parseInt(elem.data("sz"),10)))
	            elem.data("sz",elem.parent().parent().find("li.fil_val:visible").size());
	        var moreFacets = elem.next('ul').find('li');
	        moreFacets.insertBefore(elem.parent());
	    }
	    else{
	        elem.removeClass('fetch_more_exp_less');
	        elem.addClass('fetch_more_exp');
	        var sz = parseInt(elem.data("sz"),10);
	        moreFacets = elem.parent().parent().find("li.fil_val").filter(function(i){return i >= sz;});
	        elem.next().append(moreFacets);
	    }
	}
}
);
;
Portal.Portlet.Sra_Facets = Portal.Portlet.Entrez_Facets.extend ({
  
  init: function (path, name, notifier) { 
		this.base (path, name, notifier);
	},
	
	'getPortletPath': function(){
        return this.realname ;
        //return (this.realname + ".Entrez_Facets");
    }    
});

;
Portal.Portlet.Entrez_DisplayBar = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		console.info("Created DisplayBar");
		this.base(path, name, notifier);
		
		// for back button compatibility reset values when page loads
		if (this.getInput("Presentation")){
		    this.setValue("Presentation", this.getValue("LastPresentation"));
		    Portal.Portlet.Entrez_DisplayBar.Presentation = this.getValue("LastPresentation");
		}
		if (this.getInput("Format")){
		    this.setValue("Format", this.getValue("LastFormat"));
		    Portal.Portlet.Entrez_DisplayBar.Format = this.getValue("LastFormat");
		}
		if (this.getInput("PageSize")){
		    this.setValue("PageSize", this.getValue("LastPageSize"));
		    Portal.Portlet.Entrez_DisplayBar.PageSize = this.getValue("LastPageSize");
		}
		if (this.getInput("Sort")){
		    this.setValue("Sort", this.getValue("LastSort"));
		    Portal.Portlet.Entrez_DisplayBar.Sort = this.getValue("LastSort");
		}
		this.ResetDisplaySelections();
		this.ResetSendToSelection();
		
    	jQuery( 
            function(){
        
                var animationTime = jQuery("#sendto2").ncbipopper("option","openAnimationTime");
                var currentCnt = 0;
                var expTimer;
        
                function testPosition(){
                    jQuery(window).trigger("ncbipopperdocumentresize");
                    currentCnt+=10;
                    if (currentCnt<animationTime) {
                        expTimer = window.setTimeout(testPosition,10);
                    }
                }
        
                jQuery("#send_to_menu2 input").on("change click", 
                    function(){
                        currentCnt = 0;
                        if(expTimer) window.clearTimeout(expTimer);
                        testPosition();
                    } 
                );
        
            }
        );
		        
	},
	
	
	send: {
		'Cmd': null, 
		'PageSizeChanged': null,
		'ResetSendTo': null,
		'ResetCurrPage': null,
		'AddUserMessage': null
	},
		
	
	listen: {
		
		/* browser events */
			
		"sPresentation<click>": function(e, target, name){
		    this.PresentationClick(e, target, name); 
		},
		
		"sPresentation2<click>": function(e, target, name){
		    this.PresentationClick(e, target, name); 
		},
		
		"sPageSize<click>": function(e, target, name){	
		    this.PageSizeClick(e, target, name);
		},
		
		"sPageSize2<click>": function(e, target, name){	
		    this.PageSizeClick(e, target, name);
		},
		
		"sSort<click>": function(e, target, name){
		    this.SortClick(e, target, name);
		},
		
		"sSort2<click>": function(e, target, name){
		    this.SortClick(e, target, name);
		},
		
		"SetDisplay<click>": function(e, target, name){
			this.DisplayChange(e, target, name); 
		},
		
		"SendTo<click>": function(e, target, name){
			var sendto = target.value;
            var idx = target.getAttribute('sid') > 10? "2" : "";
			this.SendToClick(sendto, idx, e, target, name); 
		},
		
		"SendToSubmit<click>": function(e, target, name){
		    e.preventDefault();
		    var cmd = target.getAttribute('cmd').toLowerCase();
		    var idx = target.getAttribute('sid') > 10? "2" : "";
			this.SendToSubmitted(cmd, idx, e, target, name); 
		},
		
		/* messages from message bus*/
		
		'ResetSendTo' : function(sMessage, oData, sSrc) {
		    this.ResetSendToSelection();
		}
	
	}, // end listen
	
	
	
	/* functions */
	
	'PresentationClick': function(e, target, name){
		Portal.Portlet.Entrez_DisplayBar.Presentation = target.value;
		Portal.Portlet.Entrez_DisplayBar.Format = target.getAttribute('format');
		this.DisplayChange();
	},
	
	'PageSizeClick': function(e, target, name){ 
		Portal.Portlet.Entrez_DisplayBar.PageSize = target.value;
		this.DisplayChange();
	},
	
	'SortClick': function(e, target, name){
		Portal.Portlet.Entrez_DisplayBar.Sort = target.value;
		this.DisplayChange();
	},
	
	'DisplayChange': function(e, target, name){
	    var submit = false;
	    var extractdb = window.location.pathname.match(/\/([A-Za-z]+)\/?/); 
	    var db = (extractdb[1] && extractdb[1] != '') ? extractdb[1] : "";
	    
	    if (db != '' && getEntrezSelectedItemCount() == 1){
	        //get id, attach db and report, and link	        
	        var URL = '/' + db + '/' + getEntrezSelectedItemList() + '?report=' + Portal.Portlet.Entrez_DisplayBar.Presentation
	        + (Portal.Portlet.Entrez_DisplayBar.Format.toLowerCase() == 'text' ? '&format=text' : '');
	        window.location = URL;
	    }
	    else if (db != '' && getEntrezResultCount() == 1 && window.location.href != ""){   
	        //remove report= from URL and insert new report= into URL
	        if ((window.location.pathname != '' && window.location.pathname.match(/\/[A-Za-z]+\/\w*\d+\w*/))
	            || window.location.href.match(/\/[A-Za-z]+\/??.*term=[^&\s]+/)
	        ){
	            var URL = window.location.href.replace(/&?report=\w+/, "").replace(/\?&/, "?");
	            var hashtagindex = URL.indexOf("#");
	            if (hashtagindex >= 0){
	                URL = URL.substring(0, hashtagindex);
	            }
	            URL += (URL.match(/\?/) ? (URL.match(/\?[^\s]+/) ? "&" : "") : "?") 
	                + "report=" + Portal.Portlet.Entrez_DisplayBar.Presentation
	                + (Portal.Portlet.Entrez_DisplayBar.Format.toLowerCase() == 'text' ? '&format=text' : '');
	            window.location = URL;    
	        }
	        else {
	            submit = true;
	        }
	    }
	    else{
            submit = true;
        }
        
        if (submit){
            this.send.Cmd({'cmd': 'displaychanged'});
            
    	    this.SetPresentationChange(e, target, name);
    	    this.SetPageSizeChange(e, target, name);
    	    this.SetSortChange(e, target, name);
    	    
    	    Portal.requestSubmit();
	    }
	},
	
	'SetPresentationChange': function(e, target, name){
        this.setValue("Presentation", Portal.Portlet.Entrez_DisplayBar.Presentation);
	    this.setValue("Format", Portal.Portlet.Entrez_DisplayBar.Format);
	},
	
	'SetPageSizeChange': function(e, target, name){
	    this.setValue("PageSize", Portal.Portlet.Entrez_DisplayBar.PageSize);
		if (this.getValue("PageSize") != this.getValue("LastPageSize")){
    		//send PageSizeChanged
    		this.send.PageSizeChanged({
    			'size': this.getValue("PageSize"),
                'oldsize': this.getValue("LastPageSize")
    		});	
		}
	},
		
	'SetSortChange': function(e, target, name){
	    if (this.getInput("Sort")){
	        this.setValue("Sort", Portal.Portlet.Entrez_DisplayBar.Sort);
            if (this.getValue("Sort") != this.getValue("LastSort")){
                // ask to reset CurrPage 
    		    this.send.ResetCurrPage();
    		}
    		
    		// set sort in cookie   		
    		var extractdb = window.location.pathname.match(/\/([A-Za-z]+)\/?/); 
    	    var db = (extractdb[1] && extractdb[1] != '') ? extractdb[1] : "";
    	    
    		this.SetSortCookie(Portal.Portlet.Entrez_DisplayBar.Sort, db);
        }    	
	},
		
	'SendToClick': function(sendto, idx, e, target, name) {
		if(sendto.toLowerCase() == 'file'){
			this.SendToFile(sendto, idx);
		}
		else if(sendto.toLowerCase() == 'addtocollections'){
			this.SendToCollections(sendto, idx);
		}
		else if(sendto.toLowerCase() == 'addtoclipboard'){
		    this.SendToClipboard(sendto, idx);
		}
		else if (sendto.toLowerCase() == 'addtobibliography'){
	        this.SendToBib(sendto, e, target);
	    }
	},
	
	'SendToSubmitted': function(cmd, idx, e, target, name){
	    if (cmd == 'addtobibliography'){
	    	this.SendToBibliographySubmitted(e, cmd, idx, target);
	    }
	    else {
    	    if (cmd == 'file'){
    	         this.SendToFileSubmitted(cmd, idx, target);
    	    }
    	    else if (cmd == 'addtocollections'){
    	    	this.SendToCollectionsSubmitted(cmd, idx, target);
    	    }
    	    this.send.Cmd({'cmd': cmd});
    	    Portal.requestSubmit();
        }       
	},
	
	'ResetSendToSelection': function(){
	    var SendToInputs = this.getInputs("SendTo");
	    for (var j = 0; j < SendToInputs.length; j++){
		    if (SendToInputs[j].checked){
		        SendToInputs[j].checked = false;
			}
		}
	},
	
	'SendToFile': function(name, idx){
	    // generate content
	    var count = this.getItemCount();
		var content = 'Download ' + count + ' items.';
		this.addSendToHintContent(name, idx, content);
	},
	
	'SendToCollections': function(name, idx){
	    // generate content
        var count = this.getItemCount();
        var content= 'Add ';
        var optionNode = document.getElementById("coll_start_option" + idx);
        if (count > Portal.Portlet.Entrez_DisplayBar.CollectionsUpperLimit){
            content += Portal.Portlet.Entrez_DisplayBar.CollectionsUpperLimitText;
            if (optionNode){
            	optionNode.className = '';
            }
        }
        else{
            content += count;
            if (optionNode){
            	optionNode.className = 'hidden';
            }
        }
        content += " items.";
        this.addSendToHintContent(name, idx, content);	
	},
	
	'SendToBib': function(name, e, target){
	    jQuery('#submenu_AddToBibliography').addClass('hidden')
	    // generate content
        var count = this.getItemCount();
        var content= 'Add ';
        if (count > Portal.Portlet.Entrez_DisplayBar.BibUpperLimit){
            content += "the first " + Portal.Portlet.Entrez_DisplayBar.BibUpperLimit;
        }
        else{
            content += count;
        }
        content += " items.";
        this.addSendToHintContent(name, "", content);	
        
        // fetch other bibliography options
        var oThis = this;	
        jQuery.ui.jig.requiresLoginURL = "/account/signin/?inlinelogin=true&popuplogin=true";
        // set the menu behind myncbi login
        jQuery(target).closest('.send_to').css('z-index', '100');
        jQuery.ui.jig.requiresLogin( function(name, requiredLogin){ 
            // restore menu position after login
            jQuery(target).closest('.send_to').css('z-index', '200');
            // display message showing collection list
            jQuery('.bib_list column_list').addClass('hidden');
            jQuery('#submenu_AddToBibliography').removeClass('hidden')
            jQuery('#submenu_AddToBibliography_msg').removeClass('hidden'); 
            // fetch collection list
            var site = document.forms[0]['p$st'].value;
            xmlHttpCall(site, oThis.getPortletPath(), "GetBibliographyList", {}, oThis.ShowBibliographyList, {}, oThis);
        });
	},
	
	'SendToClipboard': function(name, idx){
	    // generate content
	    var count = this.getItemCount();
        var content= 'Add ';
        if (count > Portal.Portlet.Entrez_DisplayBar.ClipboardLimit){
            content += "the first " + Portal.Portlet.Entrez_DisplayBar.ClipboardLimit;
        }
        else{
            content += count;
        }
        content += " items.";
        this.addSendToHintContent(name, idx, content);
	},
	
	'getItemCount': function(){
	    // ask for selected items count from DbConnector
	    var selectedItemCount = getEntrezSelectedItemCount();
	    if (selectedItemCount > 0){
	        return selectedItemCount;
	    }
	    else{
	        // ask for result count from Entrez_ResultsController
	        return getEntrezResultCount();
	    }
	},
	
	'addSendToHintContent': function(name, idx, content){
	    var hintNode = document.getElementById("submenu_" + name + "_hint" + idx);
	    if (hintNode){
	        hintNode.innerHTML = content;
	        hintNode.className = 'hint';
	    }
	},
	
	'AddSendToSubmitEvent': function(){
	    // add event for SendTo submit button click. 
	    // This call is needed if the position of the submit button node has changed in relation to its parent node. 
        this.addEvent("SendToSubmit", "click", function(e, target, name) {
            var cmd = target.getAttribute('cmd');
            this.SendToSubmitted(cmd, e, target, name); 
        }, false);
    },
    
    'SendToFileSubmitted': function(cmd, idx, target){
         if (this.getInput("FFormat" + idx)){
             this.setValue("FileFormat", this.getValue("FFormat" + idx));
         }
         if (this.getInput("FSort" + idx)){
             this.setValue("FileSort", this.getValue("FSort" + idx));
         }
    },
    
    'SendToCollectionsSubmitted': function(cmd, idx, target){
         if (document.getElementById("coll_start" + idx)){
             document.getElementById("coll_startindex").value = document.getElementById("coll_start" + idx).value;
         }
    },
    
    'SendToBibliographySubmitted': function(e, cmd, idx, target){ 
        var oThis = this;	
        jQuery.ui.jig.requiresLoginURL = "/account/signin/?inlinelogin=true&popuplogin=true";
        jQuery.ui.jig.requiresLogin( function(name, requiredLogin){ 
            // update which bibliography to update
	        oThis.send.Cmd({'cmd': cmd});
	        oThis.send.AddUserMessage({'type': 'info', 
	                                    'name': 'mybib_processing_msg',
	                                    'msg': 'Adding items to bibliography ...'});
	        Portal.requestSubmit();
	        // Hack to directly make portal submit the page from async function call requiresLogin
	        // By the time asnc function finishes execution portal is not looking for submit request. It's too long after an event firing.
	        // Investigated by Mark Johnson, on 1/17/2019
	        d = Dispatcher.getInstance();
            d.submit();
	    }); 
    },
    
    'ResetDisplaySelections': function(){
        if (this.getInput("Presentation")){
            var selection = this.getValue("Presentation").toLowerCase() + this.getValue("Format").toLowerCase();
            if (document.getElementById(selection)){
                document.getElementById(selection).checked = true;
            }
            // bottom display bar
            if (document.getElementById(selection + "2")){
                document.getElementById(selection + "2").checked = true;
            }
            
        }
        if (this.getInput("PageSize")){
            var selection = 'ps' + this.getValue("PageSize");
            if (document.getElementById(selection)){
                document.getElementById(selection).checked = true;
            }
            // bottom display bar
            if (document.getElementById(selection + "2")){
                document.getElementById(selection + "2").checked = true;
            }
        }
        if (this.getInput("Sort")){
            var selection = this.getValue("Sort") || 'none'; 
            if (document.getElementById(selection)){
                document.getElementById(selection).checked = true;
            }
            // bottom display bar
            if (document.getElementById(selection + "2")){
                document.getElementById(selection + "2").checked = true;
            }
        }
    },
    
    'SetSortCookie': function(sort, db){
	    if (db != ''){
            var d = new Date();
            d.setTime(d.getTime() + (365*24*60*60*1000));
            var expires = "expires="+d.toUTCString();
            
            var newCookie = db + ":" + sort;
            var oldCookie = this.getCookie('entrezSort');
            if (oldCookie != ''){
                if (oldCookie.indexOf(db) != -1){
                    var oldSortVal = oldCookie.substring(oldCookie.indexOf(db));
                    if (oldSortVal.indexOf('&') != -1){
                        oldSortVal = oldSortVal.substring(0, oldSortVal.indexOf('&'));
                    }
                    newCookie = oldCookie.replace(oldSortVal, newCookie);
                }
                else{
                    newCookie = newCookie + "&" + oldCookie;
                }
            } 
            newCookie = "entrezSort=" + newCookie + ";domain=.ncbi.nlm.nih.gov;path=/;" + expires;
            document.cookie = newCookie;
            
		}
    },
    
    // from http://www.w3schools.com/js/js_cookies.asp
    'getCookie': function (cname) {
        var name = cname + "=";
        var ca = document.cookie.split(';');
        console.info("cookie count: " + ca.length);
        for(var i=0; i<ca.length; i++) {
            var c = ca[i];
            while (c.charAt(0)==' ') c = c.substring(1);
            if (c.indexOf(name) == 0) return c.substring(name.length,c.length);
        }
        return "";
    },
    
    'getPortletPath': function(){
        return this.realname;
    },
    
    'ShowBibliographyList': function(responseObject, userArgs){
        try {
            // Handle timeouts
            if (responseObject.status == 408) {
                //this.showMessage("Server currently unavailable. Please check connection and try again.","error");
                console.warn("Server currently unavailable. Please check connection and try again.");
            }
            else{
                var resp = '(' + responseObject.responseText + ')';
                var JSONobj = eval(resp);
                var bibList = JSONobj.BibliographyList; 
                jQuery('.mybib_list').replaceWith(bibList); 
                jQuery(".mybib_list input").on("click", function(){
                    if(jQuery(this).prop('checked')){
                        var user = jQuery(this).data('bib-id');
                        var bibName = jQuery(this).val(); 
                        jQuery('#submenu-bib-user').val(user); 
                        jQuery('#submenu-bib-name').val(bibName); 
                    }
                });
            }                       
        } catch (e) {
            //this.showMessage("Server error: " + e, "error");
            console.warn("Server error: " + e);
        }
        // in case of error, just default list will  be shown
        jQuery('#submenu_AddToBibliography_msg').addClass('hidden');  
        jQuery('.bib_list column_list').removeClass('hidden');
    }
	
},
{
    Presentation: '',
    Format: '',
    PageSize: '',
    Sort: '',
    CollectionsUpperLimit: 1000,
	CollectionsUpperLimitText: '1,000',
	ClipboardLimit: 500,
	BibUpperLimit: 200
});



;

Portal.Portlet.Sra_DisplayBar = Portal.Portlet.Entrez_DisplayBar.extend({
    
    init: function (path, name, notifier) {
        console.info("Created inherited Sra_DisplayBar");
        this.base(path, name, notifier);
        this.AddListeners(notifier);
    },
    
    // see SRA-721 (lost our ability to send Entrez search results to filter)
    SendToClick: function (sendto, idx, e, target, name) {
        if (sendto.toLowerCase() == "blast") {
            this.SelectedUids = "";
            var count = 0;
            var a = document.getElementById("ResultView");
            if (a) {
                this.SelectedUids = a.getAttribute("uid");
                count = 1;
            } else {
                var elForm = document.getElementById("EntrezForm");
                this.SelectedUids = elForm[ "EntrezSystem2.PEntrez.DbConnector.IdsFromResult"].value;
                count = getEntrezSelectedItemCount()
            }
            var content;
            var el = document.getElementById("sra_to_blast_button");
            if (count > 0) {
                content = "Send " + count + " experiment";
                if (count > 1) content += "s";
                content += ".";
                el.style.display = "";
            } else {
                content = "No results were selected";
                el.style.display = "none";
            }
            this.addSendToHintContent(sendto, idx, content);
        }
        if (sendto.toLowerCase() == "run_selector") {
            var elForm = document.getElementById("EntrezForm");
            this.SelectedUids = elForm[ "EntrezSystem2.PEntrez.DbConnector.IdsFromResult"].value;
            this.count = getEntrezSelectedItemCount();
            if (this.count > 0) {
                content = "Send " + this.count + " experiment";
                if (this.count > 1) content += "s";
                content += " to Run Selector.";
            } else {
                this.count = document.getElementById("resultcount").value;
                var el = document.getElementById("sra_run_selector_button");
                var max_count = 20000;
                if (this.count <= max_count) {
                    content = "Send whole recordset to Run Selector";
                    el.style.display = "";
                } else {
                    content = "Only first " + max_count + " experiments will be sent to Run Selector.<br/>Your recordset contains " + this.count + " experiments.";
                }
                this.count = 0; // request to run selector by webenw, not by list of UIDs
            }
            this.addSendToHintContent(sendto, idx, content);
        } else {
            this.base(sendto, idx, e, target, name);
        }
    },
    
    AddListeners: function (notifier) {
        var oThis = this;
        
        notifier.setListener(this, "sra-to-blast", function (x, y) {
            // resolve entrez uids ti sra experiments
            var oDp = new RemoteDataProvider("/entrez/eutils/esummary.fcgi?");
            oDp.onSuccess = function (obj) {
                // extract experiments
                var s = obj.responseText;
                var LEN = s.length;
                var n = "&lt;Experiment acc=&quot;";
                var len = n.length;
                var acc =[];
                while (true) {
                    var p = s.indexOf(n);
                    if (p == -1) break;
                    s = s.substring(p + len, LEN);
                    var pp = s.indexOf("&quot;");
                    if (pp == -1) break;
                    acc.push(s.substring(0, pp));
                }
                // create request to BLAST
                s = "//blast.ncbi.nlm.nih.gov/blast/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&BLAST_SPEC=SRA&DB_GROUP=Exp&NUM_ORG=";
                s += acc.length;
                p = "&EQ_MENU";
                for (var i = 0; i < acc.length;++ i) {
                    s += p + (i == 0 ? "": i) + "=" + acc[i];
                }
                var w = window.open(s);
                try {
                    w.focus();
                }
                catch (e) {
                    alert("Pop-up blocker is enabled! Please add this site to your exception list.");
                }
            };
            oDp.onError = function (obj) {
                console.info(obj.responseXML);
                alert("Error occured");
            };
            if (! oThis.SelectedUids || oThis.SelectedUids.length === 0) {
                var a = document.getElementById('maincontent').querySelectorAll('.ui-helper-hidden-accessible');
                var ids =[];
                for (var i = 0, x; x = utils.getNextSibling(a[i]);
                i++) ids.push(x.value);
                oThis.SelectedUids = ids.join(",");
            }
            oDp.Post("version=2.0&db=sra&id=" + oThis.SelectedUids);
        });
        notifier.setListener(this, "sra-to-run-selector", function (x, oData) {
            console.info(oData);
            if (oThis.count > 0) {
                window.open(oThis.run_selector_url + "uids=" + encodeURIComponent(oThis.SelectedUids));
            } else {
                window.open(oThis.run_selector_url + "WebEnv=" + oData.webenv + "&query_key=" + oData.query_key);
            }
        });
    },
    run_selector_url: "/Traces/study/?",
    count: 0,
    SelectedUids:[]
});
;
Portal.Portlet.Sra_RunSelectorPopup = Portal.Portlet.extend({
	
	init: function(path, name, notifier) {
		console.info("Created inherited Sra_RunSelectorPopup");
		this.base(path, name, notifier);	
    }
});

;
Portal.Portlet.Entrez_ResultsController = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		console.info("Created Entrez_ResultsController");
		this.base(path, name, notifier);
	},	
		
	send: {
	    'Cmd': null
	},
		
	listen: {
	
	    /* page events */
	    
	    "RemoveFromClipboard<click>": function(e, target, name){
            this.RemoveFromClipboardClick(e, target, name);
	    },
	    
		/* messages */
		
		'Cmd': function(sMessage, oData, sSrc){
		    this.ReceivedCmd(sMessage, oData, sSrc);
		},
		
		'SelectedItemCountChanged' : function(sMessage, oData, sSrc){
		    this.ItemSelectionChangedMsg(sMessage, oData, sSrc);
		},
		
		// currently sent by searchbox pubmed in journals 
		'RunLastQuery' : function(sMessage, oData, sSrc){
			if (this.getInput("RunLastQuery")){
				this.setValue ("RunLastQuery", 'true');
			}
		}
		
	},//listen
	
	'RemoveFromClipboardClick': function(e, target, name){
	    if(confirm("Are you sure you want to delete these items from the Clipboard?")){
	        this.send.Cmd({'cmd': 'deletefromclipboard'});
		    Portal.requestSubmit();  
    	}
	},
	
	// fix to not show remove selected items message when Remove from clipboard was clicked directly on one item
	'ReceivedCmd': function(sMessage, oData, sSrc){
	    if (oData.cmd == 'deletefromclipboard'){
	        Portal.Portlet.Entrez_ResultsController.RemoveOneClip = true;
	    }
	},
	
	'ItemSelectionChangedMsg': function(sMessage, oData, sSrc){
	    // do not show any messages if one item from clipbaord was removed with direct click.
	    if (Portal.Portlet.Entrez_ResultsController.RemoveOneClip){
	        Portal.Portlet.Entrez_ResultsController.RemoveOneClip = false;
	    }
	    else{
    		this.SelectedItemsMsg(oData.count);
    	    this.ClipRemoveMsg(oData.count);
    	}
	},
	
	'SelectedItemsMsg': function(count){
	    SelMsgNode = document.getElementById('result_sel');
	    if (SelMsgNode){
	        if (count > 0){
	            SelMsgNode.className = 'result_sel';
 	            SelMsgNode.innerHTML = "Selected: " + count;
 	        }
 	        else {
 	            SelMsgNode.className = 'none';
 	            SelMsgNode.innerHTML = "";
 	        }
	    }
	},
	
	'ClipRemoveMsg': function(count){
	    ClipRemNode = document.getElementById('rem_clips');
 	    if (ClipRemNode){
 	        if (count > 0){
 	            ClipRemNode.innerHTML = "Remove selected items";
 	        }
 	        else {
 	            ClipRemNode.innerHTML = "Remove all items";
 	        }
 	    }
	},
	
	'ResultCount': function(){
	    var totalCount = parseInt(this.getValue("ResultCount"));
	    totalCount = totalCount > 0 ? totalCount : 0;
	    return totalCount;
	}

},
{
    RemoveOneClip: false
});

function getEntrezResultCount() {
    var totalCount = document.getElementById("resultcount") ? parseInt(document.getElementById("resultcount").value) : 0;
	totalCount = totalCount > 0 ? totalCount : 0;
	return totalCount;
}

;
Portal.Portlet.Sra_ResultsController = Portal.Portlet.Entrez_ResultsController.extend({

	init: function(path, name, notifier) {
		this.base(path, name, notifier);
	}
	
});


;
Portal.Portlet.Entrez_Messages = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		this.base(path, name, notifier);
		
		this.setMsgAreaClassName();
	},
	
	listen: {
	   /* messages from message bus*/
		
		'AddUserMessage' : function(sMessage, oData, sSrc) {
		    // create new message node
		    var msgnode = document.createElement('li');
		    if (oData.type != ''){
		        msgnode.className = oData.type + ' icon'; 
		    }
		    if (oData.name != ''){
		        msgnode.id = oData.name; 
		    }
		    msgnode.innerHTML = "<span class='icon'>" + oData.msg + "</span>";
		    
		    // add new node as first message in message block (not ads that look like messages)
		    var parent = document.getElementById('msgportlet');
		    if (parent){
    		    var oldnode = document.getElementById(oData.name);
    		    if (oldnode){
    		        parent.removeChild(oldnode);
    		    }
    		    var firstchild = parent.firstChild;
    	        if (firstchild){
                    parent.insertBefore(msgnode, firstchild);
                }
                else{
                    parent.appendChild(msgnode);
                }
                this.setMsgAreaClassName('true');
            }
            //if there was no ul, create one, then insert the li
            else {
                var msgarea = document.getElementById('messagearea');
                if (msgarea){
                    var msgportlet = document.createElement('ul');
                    msgportlet.className = 'messages';
                    msgportlet.id = 'msgportlet';
                    msgportlet.appendChild(msgnode);
                    if (msgarea.firstChild){
                         msgarea.insertBefore(msgportlet, msgarea.firstChild);
                    }
                    else{
                        msgarea.appendChild(msgportlet);
                    }
                    this.setMsgAreaClassName('true');
                }
            }
		},
		
		'RemoveUserMessage' : function(sMessage, oData, sSrc) {
		    var msgnode = document.getElementById(oData.name);
		    if (msgnode){
		        var parent = document.getElementById('msgportlet'); 
		        if (parent){
    		        parent.removeChild(msgnode);
    		        this.setMsgAreaClassName();
    		        // if the parent ul has no children then remove the parent
    		        if (parent.firstChild){}
    		        else {
    		            if (document.getElementById('messagearea')) {
    		                document.getElementById('messagearea').removeChild(parent);
    		            }
    		        }
    		    }
		    }
		}
	}, // end listen
	
	'setMsgAreaClassName' : function(hasMsg){
        var msgarea = document.getElementById('messagearea');
	    if (msgarea){
	        var msgclass = "empty";
	        
    	    // if a message was added, hasMsg is set to true at call time to avoid checks. 
    	    // by default, hasMsg is false.
    	    if (hasMsg == 'true'){
    	        msgclass = "messagearea";
    	    }
    	    else if (msgarea.getElementsByTagName('li').length > 0){
                msgclass = "messagearea"; 
        	}
        	
            msgarea.className = msgclass;
        }
	} // end setMsgAreaClassName
});
		
		
;
Portal.Portlet.Entrez_RVBasicReport = Portal.Portlet.extend({
	
	init: function(path, name, notifier) {
		console.info("Created report portlet");
		this.base(path, name, notifier);
	},
	
	send: {
		'ItemSelectionChanged': null,
		'ClearIdList': null,
		'Cmd': null
	},
	
	listen: {
		"uid<click>" : function(e, target, name){
		    this.UidClick(e, target, name);
		},
		
		"RemoveClip<click>" : function(e, target, name){
		    this.ClipRemoveClick(e, target, name);              
		}
	},
	
	'UidClick': function(e, target, name){	
		this.send.ItemSelectionChanged( { 'id': target.value,
		                                  'selected': target.checked });
	},
	
	'ClipRemoveClick': function(e, target, name){
	    this.send.ClearIdList();
		this.send.Cmd({'cmd': 'deletefromclipboard'});
		this.send.ItemSelectionChanged( { 'id': target.getAttribute('uid'),
		                                  'selected': true });
		Portal.requestSubmit();
	}
});
   

;
Portal.Portlet.Sra_RVFull = Portal.Portlet.Entrez_RVBasicReport.extend({
	init: function(path, name, notifier) {
		console.info("Created inherited Sra_RVFull");
		this.base(path, name, notifier);
	}
});

;
(function( $ ){ // pass in $ to self exec anon fn

    // on page ready    
    
        $( 'div.portlet' ).each( function() {

            // get the elements we will need
            var $this = $( this );
            var anchor = $this.find( 'a.portlet_shutter' );
            var portBody = $this.find( 'div.portlet_content' );

            // we need an id on the body, make one if it doesn't exist already
            // then set toggles attr on anchor to point to body
            var id = portBody.attr('id') || $.ui.jig._generateId( 'portlet_content' );
            portBody.attr('id', id );
            anchor.attr('toggles', id );

            // initialize jig toggler with proper configs, then remove some classes that interfere with 
            // presentation
            var togglerOpen = anchor.hasClass('shutter_closed')? false : true; 
            anchor.ncbitoggler({
                isIcon: false,
                initOpen: togglerOpen 
            }).
                removeClass('ui-ncbitoggler-no-icon').
                removeClass('ui-widget');

            // get rid of ncbitoggler css props that interfere with portlet styling, this is hack
            // we should change how this works for next jig release
            anchor.css('position', 'absolute').
                css('padding', 0 );

            $this.find( 'div.ui-helper-reset' ).
                removeClass('ui-helper-reset');

            portBody.removeClass('ui-widget').
                css('margin', 0);

            // trigger an event with the id of the node when closed
            anchor.bind( 'ncbitogglerclose', function() {
                anchor.addClass('shutter_closed');
            });

            anchor.bind('ncbitoggleropen', function() {
                anchor.removeClass('shutter_closed');
            });

        });  // end each loop and end on page ready
})( jQuery );
/*
jQuery(document).bind('ncbitogglerclose ncbitoggleropen', function( event ) {
           var $ = jQuery;
           var eventType = event.type;
           var t = $(event.target);
           
          alert('event happened ' + t.attr('id'));
   
           if ( t.hasClass('portlet_shutter') || false ) { // if it's a portlet
               // get the toggle state
               var sectionClosed = (eventType === 'ncbitogglerclosed')? 'true' : 'false';
               alert ('now call xml-http');

            }
        });
*/

Portal.Portlet.NCBIPageSection = Portal.Portlet.extend ({
	init: function (path, name, notifier){
		this.base (path, name, notifier);
		
		this.AddListeners();
	},
    
	"AddListeners": function(){
        var oThis = this;
        
		jQuery(document).bind('ncbitogglerclose ncbitoggleropen', function( event ) {
            var $ = jQuery;
            var eventType = event.type;
            var t = $(event.target);
            
            // proceed only if this is a page section portlet {
            if ( t.hasClass('portlet_shutter')){
                var myid = '';
                if (oThis.getInput("Shutter")){
                    myid = oThis.getInput("Shutter").getAttribute('id');
                }
    
                // if the event was triggered on this portlet instance
                if (t.attr('id') && t.attr('id') == myid){
                    // get the toggle state
                    var sectionClosed = (eventType === 'ncbitogglerclose')? 'true' : 'false';
                    // react to the toggle event
                    oThis.ToggleSection(oThis.getInput("Shutter"), sectionClosed);
                }
            } // if portlet            
        });
	},
	
	"ToggleSection": function(target, sectionClosed){
	   // if remember toggle state, save the selection and log it
	   if (target.getAttribute('remembercollapsed') == 'true'){
	       this.UpdateCollapsedState(target, sectionClosed);
	   }else {
	       this.LogCollapsedState(target, sectionClosed);
	   }
	},
	
	"UpdateCollapsedState": function(target, sectionClosed){
	    var site = document.forms[0]['p$st'].value;
	    var args = { "PageSectionCollapsed": sectionClosed, "PageSectionName": target.getAttribute('pgsec_name')};
	    // Issue asynchronous call to XHR service
        var resp = xmlHttpCall(site, this.getPortletPath(), "UpdateCollapsedState", args, this.receiveCollapse, {}, this);  
	},
	
	"LogCollapsedState": function(target, sectionClosed){
	    var site = document.forms[0]['p$st'].value;
	    // Issue asynchronous call to XHR service
        var resp = xmlHttpCall(site, this.getPortletPath(), "LogCollapsedState", {"PageSectionCollapsed": sectionClosed}, this.receiveCollapse, {}, this);  
	},
	
	'getPortletPath': function(){
        return this.realname;
    }, 
    
    receiveCollapse: function(responseObject, userArgs) {
    }
	
});
		 
;
Portal.Portlet.LinkListPageSection = Portal.Portlet.NCBIPageSection.extend ({
	init: function (path, name, notifier){
		this.base (path, name, notifier);
	},
	
	"getPortletPath" : function(){
	    return (this.realname + ".NCBIPageSection");
	}
});
;
(function( $ ){ // pass in $ to self exec anon fn
    // on page ready
    $( function() {
        $('li.brieflinkpopper').each( function(){
            var $this = $( this );
            var popper = $this.find('a.brieflinkpopperctrl') ;
            var popnode = $this.find('div.brieflinkpop');
            var popid = popnode.attr('id') || $.ui.jig._generateId('brieflinkpop');
            popnode.attr('id', popid);
            popper.ncbipopper({
                destSelector: "#" + popid,
                destPosition: 'top right', 
                triggerPosition: 'middle left', 
                hasArrow: true, 
                arrowDirection: 'right',
                isTriggerElementCloseClick: false,
                adjustFit: 'none',
                openAnimation: 'none',
                closeAnimation: 'none',
                delayTimeout : 130
            });
        }); // end each loop  
    });// end on page ready
})( jQuery );

Portal.Portlet.BriefLinkPageSection = Portal.Portlet.LinkListPageSection.extend({

	init: function(path, name, notifier) {
	    console.info("Created BriefLinkPageSection");
		this.base(path, name, notifier);
	},
	
	"getPortletPath" : function(){
	    return (this.realname + ".LinkListPageSection.NCBIPageSection");
	}
	
});
;
Portal.Portlet.DiscoveryDbLinks = Portal.Portlet.BriefLinkPageSection.extend({
    
    init: function(path, name, notifier) {
		this.base(path, name, notifier);
	},
	
	"getPortletPath" : function(){
	    return (this.realname + ".BriefLinkPageSection.LinkListPageSection.NCBIPageSection");
	}
});
;
(function( $ ){ // pass in $ to self exec anon fn
    // on page ready
    $( function() {
        $('li.ralinkpopper').each( function(){
            var $this = $( this );
            var popper = $this;
            var popnode = $this.find('div.ralinkpop');
            var popid = popnode.attr('id') || $.ui.jig._generateId('ralinkpop');
            popnode.attr('id', popid);
            popper.ncbipopper({
                destSelector: "#" + popid,
                destPosition: 'top right', 
                triggerPosition: 'middle left', 
                hasArrow: true, 
                arrowDirection: 'right',
                isTriggerElementCloseClick: false,
                adjustFit: 'none',
                openAnimation: 'none',
                closeAnimation: 'none',
                delayTimeout : 130
            });
        }); // end each loop  
    });// end on page ready
})( jQuery );

Portal.Portlet.HistoryDisplay = Portal.Portlet.NCBIPageSection.extend({

	init: function(path, name, notifier) {
		console.info("Created History Ad...");
		this.base(path, name, notifier);    
	},
	
	send: {
      'Cmd': null      
    },   
    
    receive: function(responseObject, userArgs) {  
         var cmd = userArgs.cmd;
         var rootNode = document.getElementById('HTDisplay'); 
         var ul = document.getElementById('activity');
         var resp = responseObject.responseText;
             
         if (cmd == 'HTOn') { 
            rootNode.className = '';    // hide all msg and the turnOn link
            try {
            //alert(resp);
                // Handle timeouts
                if (responseObject.status == 408) { 
                    rootNode.className = 'HTOn'; // so that the following msg will show up
                    rootNode.innerHTML = "<p class='HTOn'>Your browsing activity is temporarily unavailable.</p>";
                    return;
                }
                   
                 // Looks like we got something...
                 resp = '(' + resp + ')';
                 var JSONobj = eval(resp);
                 
                 // Build new content (ul)
                 var newHTML = JSONobj.Activity;
                 var newContent = document.createElement('div');
                 newContent.innerHTML = newHTML;
                 var newUL = newContent.getElementsByTagName('ul')[0];
                 //alert(newHTML);
                 //alert(newContent.innerHTML);
                 //alert(newUL.innerHTML);
                 // Update content
                 rootNode.replaceChild(newUL, ul);
                 //XHR returns no activity (empty ul), e.g. activity cleared
                 if (newUL.className == 'hide')                     
                     rootNode.className = 'HTOn';  // show "Your browsing activity is empty." message
                 
            }         
            catch (e) {
                //alert('error');
                rootNode.className = 'HTOn'; // so that the following msg will show up
                rootNode.innerHTML = "<p class='HTOn'>Your browsing activity is temporarily unavailable.</p>";
           }
         }
         else if (cmd == 'HTOff') {                         
             if (ul != null) { 
                 ul.className='hide'; 
                 ul.innerHTML = ''; // clear activity
             }
             rootNode.className = 'HTOff';    // make "Activity recording is turned off." and the turnOn link show up             
         }
         else if (cmd == 'ClearHT') { 
             var goAhead = confirm('Are you sure you want to delete all your saved Recent Activity?');
             if (goAhead == true) { 
                 if ( rootNode.className == '') { //                 
                     rootNode.className = 'HTOn';  // show "Your browsing activity is empty." message                                  
                     if (ul != null) {
                         ul.className='hide'; 
                         ul.innerHTML = '';
                     }
                 }
             }
         } 
         
    },
    
	listen: {
	  'Cmd' : function(sMessage, oData, sSrc){
			console.info("Inside Cmd in HistoryDisplay: " + oData.cmd);
			this.setValue("Cmd", oData.cmd);
	  },	  
		
      "HistoryToggle<click>" : function(e, target, name){
         //alert(target.getAttribute("cmd"));
         this.send.Cmd({'cmd': target.getAttribute("cmd")});         
         console.info("Inside HistoryToggle in HistoryDisplay: " + target.getAttribute("cmd"));
         
         //var site = document.forms[0]['p$st'].value;
         var cmd =  target.getAttribute("cmd");     
               
         // Issue asynchronous call to XHR service, callback is to update the portlet output            
         this.doRemoteAction(target.getAttribute("cmd"));                      
      }, 
      
      "HistoryOn<click>" : function(e, target, name){
         this.send.Cmd({'cmd': target.getAttribute("cmd")});
         //$PN('Pubmed_ResultsSearchController').getInput('RecordingHistory').value = 'yes';		 
         console.info("Inside HistoryOn in HistoryDisplay: " + target.getAttribute("cmd"));
         this.doRemoteAction(target.getAttribute("cmd"));         
      },
      
      "ClearHistory<click>" : function(e, target, name){
         this.send.Cmd({'cmd': target.getAttribute("cmd")});
         this.doRemoteAction(target.getAttribute("cmd"));         
      }
    },
    
    'getPortletPath': function(){
        return this.realname + ".NCBIPageSection";
    }, 
    
    'doRemoteAction': function(command) {
         var site = document.forms[0]['p$st'].value;          
	     var resp = xmlHttpCall(site, this.realname, command, {}, this.receive, {'cmd': command}, this);
    }
});

;
Portal.Portlet.DbConnector = Portal.Portlet.extend({

	init: function(path, name, notifier) {
		var oThis = this;
		console.info("Created DbConnector");
		this.base(path, name, notifier);
		
		// reset Db value to original value on page load. Since LastDb is the same value as Db on page load and LastDb is not changed on
		// the client, this value can be used to reset Db. This is a fix for back button use.
		if (this.getValue("Db") != this.getValue("LastDb")){
		    this.setValue("Db", this.getValue("LastDb"));
		}
     
		// the SelectedIdList and id count from previous iteration (use a different attribute from IdsFromResult to prevent back button issues)
		Portal.Portlet.DbConnector.originalIdList = this.getValue("LastIdsFromResult");
		console.info("originalIdList " + Portal.Portlet.DbConnector.originalIdList);
		// if there is an IdList from last iteration set the count
		if (Portal.Portlet.DbConnector.originalIdList != ''){
			Portal.Portlet.DbConnector.originalCount = Portal.Portlet.DbConnector.originalIdList.split(/,/).length;
		}

		notifier.setListener(this, 'HistoryCmd', 
        	function(oListener, custom_data, sMessage, oNotifierObj) {
           		var sbTabCmd = $N(oThis.path + '.TabCmd');
           		sbTabCmd[0].value = custom_data.tab;
        	}
    		, null);
    
	},

	send: {
   		'SelectedItemCountChanged': null,
   		'newUidSelectionList': null,
   		'SavedSelectedItemCount': null,
   		'SavedUidList': null
	},

	listen: {
	
		//message from Display bar on Presentation change 
		'PresentationChange' : function(sMessage, oData, sSrc){
			
			// set link information only if it exists
			if (oData.dbfrom){
				console.info("Inside PresentationChange in DbConnector: " + oData.readablename);
				this.setValue("Db", oData.dbto);
				this.setValue("LinkSrcDb", oData.dbfrom);
				this.setValue("LinkName", oData.linkname);
				this.setValue("LinkReadableName", oData.readablename);
			}
			//document.forms[0].submit();
		},
		
		// various commands associated with clicking different form control elements
		'Cmd' : function(sMessage, oData, sSrc){
			console.info("Inside Cmd in DbConnector: " + oData.cmd);
			this.setValue("Cmd", oData.cmd);
			
			// back button fix, clear TabCmd
			if (oData.cmd == 'Go' || oData.cmd == 'PageChanged' || oData.cmd == 'FilterChanged' || 
			oData.cmd == 'DisplayChanged' || oData.cmd == 'HistorySearch' || oData.cmd == 'Text' || 
			oData.cmd == 'File' || oData.cmd == 'Printer' || oData.cmd == 'Order' || 
			oData.cmd == 'Add to Clipboard' || oData.cmd == 'Remove from Clipboard' || 
			oData.cmd.toLowerCase().match('details')){
				this.setValue("TabCmd", '');
				console.info("Inside Cmd in DbConnector, reset TabCmd: " + this.getValue('TabCmd'));
			}

		},
		
		
		// the term to be shown in the search bar, and used from searching
		'Term' : function(sMessage, oData, sSrc){
			console.info("Inside Term in DbConnector: " + oData.term);
			this.setValue("Term", oData.term);
		},
		
		
		// to indicate the Command Tab to be in
		'TabCmd' : function(sMessage, oData, sSrc){
			console.info("Inside TABCMD in DbConnector: " + oData.tab);
			this.setValue("TabCmd", oData.tab);
			console.info("DbConnector TabCmd: " + this.getValue("TabCmd"));
		},
		
		
		// message sent from SearchBar when db is changed while in a Command Tab
		'DbChanged' : function(sMessage, oData, sSrc){
			console.info("Inside DbChanged in DbConnector");
			this.setValue("Db", oData.db);
		},
		
		// Handles item select/deselect events
		// Argument is { 'id': item-id, 'selected': true or false }
		'ItemSelectionChanged' : function(sMessage, oData, oSrc) {
			var sSelection = this.getValue("IdsFromResult");
			var bAlreadySelected = (new RegExp("\\b" + oData.id + "\\b").exec(sSelection) != null);
	       	var count =0;
	       	
			if (oData.selected && !bAlreadySelected) {
				sSelection += ((sSelection > "") ? "," : "") + oData.id;
			   	this.setValue("IdsFromResult", sSelection);
			   	if (sSelection.length > 0){
			   		count = sSelection.split(',').length;
			   	}
			   	this.send.SelectedItemCountChanged({'count': count});
			   	this.send.newUidSelectionList({'list': sSelection});
			   	jQuery(document).trigger("itemsel",{'list': sSelection});
		   	} else if (!oData.selected && bAlreadySelected) {
				sSelection = sSelection.replace(new RegExp("^"+oData.id+"\\b,?|,?\\b"+oData.id+"\\b"), '');
		   	   	this.setValue("IdsFromResult", sSelection);
				console.info("Message ItemSelectionChanged - IdsFromResult after change:  " + this.getValue("IdsFromResult"));
			   	if (sSelection.length > 0){
			   		count = sSelection.split(',').length;
			   	}
				console.info("Message ItemSelectionChanged - IdsFromResult length:  " + count);   
				this.send.SelectedItemCountChanged({'count': count});
			   	this.send.newUidSelectionList({'list': sSelection});
			   	jQuery(document).trigger("itemsel",{'list': sSelection});
		   	}
		},
				
		// FIXME: This is the "old message" that is being phased out.
		// when result citations are selected, the list of selected ids are intercepted here,
		// and notification sent that selected item count has changed.
		'newSelection' : function(sMessage, oData, sSrc){
		
			// Check if we already have such IDs in the list
			var newList = new Array();
			var haveNow = new Array();
			if(Portal.Portlet.DbConnector.originalIdList){
				haveNow = Portal.Portlet.DbConnector.originalIdList.split(',');
				newList = haveNow;
			}
			
			var cameNew = new Array();
			if (oData.selectionList.length > 0) {
				cameNew = oData.selectionList;
			}
			
			if (cameNew.length > 0) {
				for(var ind=0;ind<cameNew.length;ind++) {
					var found = 0;
					for(var i=0;i<haveNow.length;i++) {
						if (cameNew[ind] == haveNow[i]) {
							found = 1;
							break;
						}
					}
						//Add this ID if it is not in the list
					if (found == 0) {
						newList.push(cameNew[ind]);
					}
				}
			}
			else {
				newList = haveNow;
			}

				// if there was an IdList from last iteration add new values to old
			var count = 0;
			if ((newList.length > 0) && (newList[0].length > 0)){
				count = newList.length;
			}
			
			console.info("id count = " + count);
			this.setValue("IdsFromResult", newList.join(","));
			
			this.send.SelectedItemCountChanged({'count': count});
			this.send.newUidSelectionList({'list': newList.join(",")});
			jQuery(document).trigger("itemsel",{'list': newList.join(",")});
		},


		// empty local idlist when list was being collected for other purposes.
		//used by Mesh and Journals (empty UidList should not be distributed, otherwise Journals breaks)
		// now used by all reports for remove from clipboard function.
		'ClearIdList' : function(sMessage, oData, sSrc){
			this.setValue("IdsFromResult", '');
			this.send.SelectedItemCountChanged({'count': '0'});
			this.send.newUidSelectionList({'list': ''});
			jQuery(document).trigger("itemsel",{'list': ""});
		}, 


		// back button fix: when search backend click go or hot enter on term field,
		//it also sends db. this db should be same as dbconnector's db
		'SearchBarSearch' : function(sMessage, oData, sSrc){
			if (this.getValue("Db") != oData.db){
				this.setValue("Db", oData.db);
			}
		},
		
		// back button fix: whrn links is selected from DisplayBar,
		//ResultsSearchController sends the LastQueryKey from the results on the page
		// (should not be needed by Entrez 3 code)
		'LastQueryKey' : function(sMessage, oData, sSrc){
			if (this.getInput("LastQueryKey")){
				this.setValue("LastQueryKey", oData.qk);
			}
		},
		
		'QueryKey' : function(sMessage, oData, sSrc){
			if (this.getInput("QueryKey")){
				this.setValue("QueryKey", oData.qk);
			}
		},
		
		
		//ResultsSearchController asks for the initial item count in case of send to file 
		'needSavedSelectedItemCount' : function(sMessage, oData, sSrc){
			var count = 0;
			if(this.getInput("IdsFromResult")){
				if (this.getValue("IdsFromResult").length > 0){
					count = this.getValue("IdsFromResult").split(',').length;
				}
				console.info("sending SavedSelectedItemCount from IdsFromResult: " + count);
			}
			else{
				count = Portal.Portlet.DbConnector.originalCount;
				console.info("sending SavedSelectedItemCount from OriginalCount: " + count);
			}
			this.send.SavedSelectedItemCount({'count': count});
		},
		
		// Force form submit, optionally passing db, term and cmd parameters
		'ForceSubmit': function (sMessage, oData, sSrc)
		{
		    if (oData.db)
    			this.setValue("Db", oData.db);
		    if (oData.cmd)
    			this.setValue("Cmd", oData.cmd);
		    if (oData.term)
    			this.setValue("Term", oData.term);
    		Portal.requestSubmit ();
		},
		
		'LinkName': function (sMessage, oData, sSrc){
		    this.setValue("LinkName", oData.linkname);
		},
		
		'IdsFromResult': function (sMessage, oData, sSrc){
		    this.setValue("IdsFromResult", oData.IdsFromResult);
		},
		
		'SendSavedUidList': function (sMessage, oData, sSrc){
		    this.send.SavedUidList({'idlist': this.getValue("IdsFromResult")});
		}
		
	}, //listen
	
	/* other portlet functions */
	
	// DisplayBar in new design wants selected item count
	'SelectedItemCount': function(){
	    var count = 0;
		if(this.getInput("IdsFromResult")){
			if (this.getValue("IdsFromResult") != ''){
				count = this.getValue("IdsFromResult").split(',').length;
			}
		}
		else{
			count = Portal.Portlet.DbConnector.originalCount;
		}
		return count;
	},
	
	'SelectedItemList': function(){
		if(this.getInput("IdsFromResult") && this.getValue("IdsFromResult") != ''){
			return this.getValue("IdsFromResult");
		}
		else{
			return Portal.Portlet.DbConnector.originalIdList;
		}
		
	},
	setValue: function(name, value){
	    if(name == 'Term')
	        value = jQuery.trim(value);
	    this.base(name,value);
	}
},
{
	originalIdList: '',
	originalCount: 0
});

function getEntrezSelectedItemCount() {
    return $PN('DbConnector').SelectedItemCount();
}

function getEntrezSelectedItemList() {
    return $PN('DbConnector').SelectedItemList();
}
