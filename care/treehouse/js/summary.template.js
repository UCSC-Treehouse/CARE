// populateSections() is run when the page is loaded.
// It calls populate() on normal items and parseCustom() on doCustomParse items.
// parseCustom calls other functions on an item-specific basis
// requires an 'annotations' data structure with items in the format 'elementID' : { 'value' : ... }

    // Populate the document with the analyst-edited text in the 'annotations' data structure
      function populateSections() {
          
        for(let item in annotations){
          // If an item has a truthy doCustomParse key, parseCustom instead
          if (annotations[item]["doCustomParse"]){
              parseCustom(item, annotations[item])
          } else {
              populate(item, annotations[item]["value"]);
          }
        }
        // Then, resize the slides' imgcaption items based on how much (post-population) text is in them.
        // Algorithm was calculated empirically based on the following observations and is at minimum 12 pt.
        // 100 chars: 35pt; 200 chars: 24pt; 300 chars: 18pt; 400 chars; 16pt
        if(documentTypeIsSlides){
          for(let elem of document.getElementsByClassName("imgcaption")){
              
              let charcount = elem.innerHTML.length;
              let calcfont = Math.round(2**(5.5 - 0.01 * charcount) + 12)
              elem.style.fontSize = calcfont + "px";
          }
        }
      } // end populateSections
    
    // Helper - If an element is present in the document, set innerHTML to text
    function populate(elementId, text){
        if(document.getElementById(elementId) !== null){
            document.getElementById(elementId).innerHTML=text;
        } else {
            console.log("Couldn't find '" + elementId + "' in the document to write contents.");
        }
    }


    // Parse fields that are more than simply populating with the value 
    // The custom parser function that is named "item" will be called, and the 'contents' dict passed in.
    function parseCustom(item, contents){
        custom_parsers = {
            "hideable_items" : function(info){
                hideSelectedHideables(info["value"]);
            },
            "molecular_category" : function(info){
                parseMolecularCategories(info["value"], info["value_of_other"]);
            },
            "mcs_mutations_table" : function(info){
                updateMutationsTable(info["value"]);
            },
            "mcs_mutations_table_diagnosis" : function(info){
                updateMutationsTable_diagnosis(info["value"]);
            },
            "prior_molecular_testing" : function(info){
                parsed_item = addPriorMolecularTesting(info["value"]);
                populate("prior_molecular_testing", parsed_item);
            },
            "analyst_custom_images" : function(info){
                parsed_item = addCustomImageSlides(info["value"]);
                populate("analyst_custom_images", parsed_item);
            },
            "all_molecular_findings" : function(info){
                parsed_item = parseAllMolecularFindings(info["value"]);
                populate("all_molecular_findings", parsed_item);
                
                // Also populate standalone slides or summary references list
                if(documentTypeIsSlides){
                    slides_text = parseAllMolecularFindings_slides(info["value"]);
                    populate("all_findings_individual_slides", slides_text);
                } else {
                    summary_text = allMolecularFindingsReferences_summary(info["value"]);
                    populate("finding_specific_references", summary_text);
               }
            },
            "custom_font_size_personalized_cohorts_slide" : function(info){
                size = info["value"] // string - so empty or otherwise
                if(documentTypeIsSlides && size){
                    document.getElementById("tablePersonalizedCohortCounts").style.fontSize = size;
                }
            },
            "quality_control_statement": function (info){
                fail_reason = info["value"] // string - so empty or otherwise
                if((!documentTypeIsSlides) && fail_reason){
                    populate("quality_control_statement", fail_reason);
                    // Also make the box red
                    document.getElementById("highlight_failed_qc_div").classList.add("red");
                }
            }
        } // end custom_parsers dict
        custom_parsers[item](contents); // Call the parser
    }
 

      // Some items should only be displayed if this is a Stanford Registry sample
      // Usage: add class registryOnly or noRegistry to things.
      function toggleRegistryItems(){
          let isRegistry = (annotations["generic_or_registry"]["value"] == "registry");
          let showHide = ["registryOnly", "noRegistry"] // Show first, hide second
          
          if(! isRegistry) {
              showHide = showHide.reverse();
          }
          
          for(let elem of document.getElementsByClassName(showHide[0])){
              elem.classList.remove("hiddenForRegistryToggle");
          }
          for(let elem of document.getElementsByClassName(showHide[1])){
              elem.classList.add("hiddenForRegistryToggle");
          }
      }

      // Takes array of object IDs to hide
      function hideSelectedHideables(hideables){
          for (let hideMe of hideables){
              if(document.getElementById(hideMe)){
                  document.getElementById(hideMe).classList.add("displayNone");
              }
          }
      }
    
      function addPriorMolecularTesting(rows){
          table = "";
          
          for (let row of rows){
              table += "<tr>";
              table += '<td style="padding:0px;">';
              table += row;
              table +='</td>';
              table +='</tr>';
          }
          return table;
      }
      // mutation key = "mutations_SAMPLEID", value = string comma-separated desired mutations.
      // Will only replace if nonempty; to clear, leave it as a space.
      function updateMutationsTable(mutations){
          for (let mutation in mutations){
              if(mutations[mutation].length > 0){
                  if(document.getElementById(mutation)){
                     document.getElementById(mutation).innerHTML = mutations[mutation];
                     }
              }
          }
      }
      // diagnosis key = "diagnosis_SAMPLEID", value = string diagnosis.
      // Will only replace if nonempty; to clear, leave it as a space.
      function updateMutationsTable_diagnosis(diagnoses){
          for (let diagnosis in diagnoses){
              if(diagnoses[diagnosis].length > 0){
                  if(document.getElementById(diagnosis)){
                     document.getElementById(diagnosis).innerHTML = diagnoses[diagnosis];
                     }
              }
          }
      }

      function allMolecularFindingsReferences_summary(rows) {
          // Get all the references
          // Sort by key
          // Populate the OL
          all_references = {}
          for (let row of rows) {
              all_references = Object.assign(all_references, row["references"])
          }
          reference_numbers = Object.keys(all_references);
          reference_numbers.sort(function(a, b){return a-b});
          list_items = ""
          for (let ref of reference_numbers){
              list_items += '<li>';
              list_items += all_references[ref];
              list_items += "</li>"
          }
          return list_items;
      }

      // Make a page div for each extra image desired
      function addCustomImageSlides(rows){
          let pages = "";
          for (let row of rows){
              // If location is tumormap or pathway, replace that slide + summary item instead of adding a new one
              if( (row['location'] == 'tumormap') || (row['location'] == 'pathway')){
                  let boxname = row['location'];
                  let title = document.getElementById(boxname + "_title");
                  let img = document.getElementById(boxname + "_img");
                  let caption = document.getElementById(boxname + "_caption");
                  
                  title.innerHTML = row['title'];
                  img.src = row['img'];
                  caption.innerHTML = row['caption'];
                  
              } else {
                  pages += '<div class="page">';
                  pages += '<h4 class="ui top center aligned inverted header">' + row['title'] + '</h4>';
                  pages += '<div><img src="' + row["img"] + '" style="width:9.9in;height:5.5in;">';
                  pages += '<div class="ui center very basic segment imgcaption">'
                  pages += row["caption"];
                  pages +=  '</div></div></div>';
              }
          }
          return pages;
      }
        
      // Make a page div for each finding
      function parseAllMolecularFindings_slides(rows){
          let pages = "";
          let legend_text = document.getElementById("findingsLegend").innerHTML;
          // Can't use column-count so split the legend_text in two
          legend_text = '<div>' + legend_text + '</div>';
          legend_text = legend_text.replace('<li value="e)"', 
                                            '</ol></div><div style="padding-left:1%;"><ol class="ui list" ><li value="e)"');
          

          for (let row of rows){
              if(row['hide'] == 'yes'){
                  continue; // Don't add markup for hidden rows
              }
              else {
                  pages += '<div class="page">';
                  pages += '<h4 class="ui top center aligned inverted header">Molecular Findings</h4>';
                  
                  if( Object.keys(row["references"]).length > 0){
                      pages += '<div class="ui top attached segment" style="margin:0px;">';
                  } else {
                      pages += '<div class="ui segment" style="margin:0px;">';
                  }
                  pages += '<span class="ui horizontal label" style="margin-top:10px;">Molecular Abnormality</span>' 
                      + row["abnormality"] + "<br/>" ;
                  pages += '<span class="ui horizontal label" style="margin-top:10px;">Category</span>' 
                      + row["category"] + "<sup>" + row["supports"] + "</sup>" + "<br/>" ;
                  pages += '<span class="ui horizontal label" style="margin-top:10px;">Associated Drugs</span>' 
                      + row["therapy"] + "<br/>" ;
                  pages +=  row["analyst_summary"];
   
                  if( Object.keys(row["references"]).length > 0)
                  {
                      pages += '</div><div class="ui bottom attached segment"';
                      pages += 'style="font-size:60%;line-height:1.2;margin-bottom:1%;">';
                      pages += '<span class="ui horizontal label" style="margin-top:10px;">References</span><br/>' ;
                      for(ref in row["references"]){
                      pages += row["references"][ref] + "<br/>";
                      }
                  }
                  pages += '</div>';
                  pages += '<div style="margin:0px;padding:5px;font-size:50%;display:flex;zoom:80%">';
                  pages += legend_text;
                  pages += '</div>';

                  pages += '</div>'; // end segment div, and page div
              }
          }
          return pages;
      }

      // Takes a list of row objects and htmlifies it
      // For slides, skip the analyst summary
      function parseAllMolecularFindings(rows){
          table = "";
          for (let row of rows){
              table += "<tr>";
              table += "<td>" + row["abnormality"] + "</td>";
              table += "<td>" + row["category"]+ "<sup>" + row["supports"] + "</sup>" + "</td>";
              table += "<td>" + row["therapy"] + "</td>";
              if(! documentTypeIsSlides){
              table += "<td>" + row["analyst_summary"] + "</td>";
              }
              table += "</tr>";
          }
          return table;
      }
      // Bold the chosen categories and also check their boxes. Also fill in "other" if present
      function parseMolecularCategories(categories, other_value){
          for(let category in categories){
              if(categories[category]){
               document.getElementById(category).classList.add("chosen_molecular_category");
               document.getElementById(category).childNodes[0].checked = true;
               if(category == "other"){
                  document.getElementById(category).childNodes[1].nodeValue = " " + other_value;
               }
              }
          }
      }
      // Today as MM/DD/YYYY , offset for current timezone
      function todaysDate(){
        now = new Date();
        t=new Date(now.getTime() - (now.getTimezoneOffset() * 60000)).toJSON().slice(0,10).split("-");
        today=[t[1],t[2],t[0]].join("/");
        return today;
      }
