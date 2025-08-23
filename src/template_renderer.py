"""
Template rendering utilities for HTML report generation.
"""
import os
from typing import Dict, List
import html

class HTMLTemplateRenderer:
    """Handles HTML template rendering for Shiba reports."""
    
    def __init__(self, template_dir: str):
        self.template_dir = template_dir
        
    def load_template(self, filename: str) -> str:
        """Load a template file."""
        template_path = os.path.join(self.template_dir, filename)
        with open(template_path, 'r', encoding='utf-8') as f:
            return f.read()
    
    def escape_and_format_iframe_content(self, lines: List[str]) -> str:
        """Format iframe content from plotly HTML lines."""
        return '\n'.join([html.escape(line.strip()) for line in lines])
    
    def render_splicing_section(self, event_data: Dict) -> str:
        """Render a single splicing event section."""
        template = self.load_template('splicing_section.html')
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{section_id}', str(event_data['id']))
        result = result.replace('{icon_class}', str(event_data['icon']))
        result = result.replace('{section_title}', str(event_data['title']))
        result = result.replace('{section_description}', str(event_data['description']))
        result = result.replace('{volcano_content}', str(event_data['volcano_content']))
        result = result.replace('{scatter_content}', str(event_data['scatter_content']))
        
        return result
    
    def render_summary_html(self, data: Dict) -> str:
        """Render complete summary HTML (requires external CSS/JS)."""
        template = self.load_template('summary.html')
        
        # Generate all splicing event sections
        sections = []
        for event in data['splicing_events']:
            section_html = self.render_splicing_section(event)
            sections.append(section_html)
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{shiba_command}', str(data['shiba_command']))
        result = result.replace('{pca_tpm_content}', str(data['pca_tpm_content']))
        result = result.replace('{pca_psi_content}', str(data['pca_psi_content']))
        result = result.replace('{splicing_sections}', '\n'.join(sections))
        
        return result

    def render_standalone_summary_html(self, data: Dict) -> str:
        """Render standalone summary HTML with embedded CSS/JS."""
        template = self.load_template('standalone_summary.html')
        
        # Generate all splicing event sections (standalone versions)
        sections = []
        for event in data['splicing_events']:
            section_html = self.render_splicing_section_standalone(event)
            sections.append(section_html)
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{shiba_command}', str(data['shiba_command']))
        result = result.replace('{pca_tpm_content}', str(data['pca_tpm_content']))
        result = result.replace('{pca_psi_content}', str(data['pca_psi_content']))
        result = result.replace('{splicing_sections}', '\n'.join(sections))
        
        return result

    def render_splicing_section_standalone(self, event_data: Dict) -> str:
        """Render a single splicing event section for standalone HTML."""
        # For standalone version, we embed everything inline
        return f'''
            <section class="splicing-section" id="{event_data['id']}">
                <div class="section-header">
                    <div class="section-icon">
                        <i class="{event_data['icon']}"></i>
                    </div>
                    <div class="section-info">
                        <h2 class="section-title">{event_data['title']}</h2>
                        <p class="section-description">{event_data['description']}</p>
                    </div>
                </div>
                
                <div class="plots-grid">
                    <div class="plot-container">
                        <h3 class="plot-title">
                            <i class="fas fa-volcano"></i>
                            Volcano Plot
                        </h3>
                        <div class="plot-frame">
                            {event_data['volcano_content']}
                        </div>
                        <p class="plot-description">Differential splicing significance vs. magnitude</p>
                    </div>
                    
                    <div class="plot-container">
                        <h3 class="plot-title">
                            <i class="fas fa-chart-scatter"></i>
                            Scatter Plot
                        </h3>
                        <div class="plot-frame">
                            {event_data['scatter_content']}
                        </div>
                        <p class="plot-description">PSI correlation between conditions</p>
                    </div>
                </div>
            </section>
        '''

    def render_individual_event_html(self, event_data: Dict) -> str:
        """Render HTML for an individual splicing event (requires external CSS/JS)."""
        template = self.load_template('individual_event.html')
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{event_title}', str(event_data['title']))
        result = result.replace('{event_icon}', str(event_data['icon']))
        result = result.replace('{event_description}', str(event_data['description']))
        result = result.replace('{volcano_content}', str(event_data['volcano_content']))
        result = result.replace('{scatter_content}', str(event_data['scatter_content']))
        
        return result

    def render_event_card(self, event_data: Dict) -> str:
        """Render an event card for the index page."""
        return f'''
            <a href="data/{event_data['id']}.html" class="event-card">
                <div class="event-header">
                    <div class="event-icon">
                        <i class="{event_data['icon']}"></i>
                    </div>
                    <h3 class="event-title">{event_data['title']}</h3>
                </div>
                <p class="event-description">{event_data['description']}</p>
                <div class="event-stats">
                    <span class="view-link">
                        View Details
                        <i class="fas fa-arrow-right"></i>
                    </span>
                    <span class="event-count">
                        {event_data.get('event_count', 0)} events
                    </span>
                </div>
            </a>
        '''
    
    def render_index_html(self, data: Dict) -> str:
        """Render the main index/overview HTML (requires external CSS/JS)."""
        template = self.load_template('index.html')
        
        # Render all event cards
        event_cards = []
        for event in data['splicing_events']:
            card_html = self.render_event_card(event)
            event_cards.append(card_html)
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{shiba_command}', str(data['shiba_command']))
        result = result.replace('{shiba_version}', str(data.get('shiba_version', 'unknown')))
        result = result.replace('{pca_tpm_content}', str(data['pca_tpm_content']))
        result = result.replace('{pca_psi_content}', str(data['pca_psi_content']))
        result = result.replace('{splicing_summary_content}', str(data.get('splicing_summary_content', '')))
        result = result.replace('{event_cards}', '\n'.join(event_cards))
        
        return result

    def render_index_standalone_html(self, data: Dict) -> str:
        """Render standalone index HTML with embedded CSS/JS."""
        template = self.load_template('index_standalone.html')
        
        # Generate event cards HTML
        event_cards = []
        for event in data['splicing_events']:
            card_html = self.render_event_card(event)
            event_cards.append(card_html)
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{shiba_command}', str(data['shiba_command']))
        result = result.replace('{shiba_version}', str(data.get('shiba_version', 'unknown')))
        result = result.replace('{pca_tpm_content}', str(data['pca_tpm_content']))
        result = result.replace('{pca_psi_content}', str(data['pca_psi_content']))
        result = result.replace('{splicing_summary_content}', str(data.get('splicing_summary_content', '')))
        result = result.replace('{event_cards}', '\n'.join(event_cards))
        
        return result

    def render_individual_event_standalone_html(self, event_data: Dict) -> str:
        """Render standalone HTML for an individual splicing event (completely self-contained)."""
        template = self.load_template('individual_event_standalone.html')
        
        # Use replace instead of format to avoid conflicts with CSS braces
        result = template
        result = result.replace('{event_title}', str(event_data['title']))
        result = result.replace('{event_icon}', str(event_data['icon']))
        result = result.replace('{event_description}', str(event_data['description']))
        result = result.replace('{volcano_content}', str(event_data['volcano_content']))
        result = result.replace('{scatter_content}', str(event_data['scatter_content']))
        
        return result

    def copy_static_files(self, output_dir: str):
        """Copy CSS and JS files to output directory."""
        import shutil
        
        static_files = ['styles.css', 'script.js']
        for filename in static_files:
            src = os.path.join(self.template_dir, filename)
            dst = os.path.join(output_dir, filename)
            if os.path.exists(src):
                shutil.copy2(src, dst)

def get_splicing_event_config():
    """Return configuration for splicing event types."""
    return [
        {
            'id': 'se',
            'code': 'SE',
            'title': 'Skipped Exon (SE)',
            'icon': 'fas fa-minus-square',
            'description': 'Events where an exon is either included or skipped.'
        },
        {
            'id': 'five',
            'code': 'FIVE',
            'title': 'Alternative 5′ Splice Site (FIVE)',
            'icon': 'fas fa-arrow-left',
            'description': 'Events involving alternative donor splice sites at the 5′ end of introns.'
        },
        {
            'id': 'three',
            'code': 'THREE',
            'title': 'Alternative 3′ Splice Site (THREE)',
            'icon': 'fas fa-arrow-right',
            'description': 'Events involving alternative acceptor splice sites at the 3′ end of introns.'
        },
        {
            'id': 'mxe',
            'code': 'MXE',
            'title': 'Mutually Exclusive Exons (MXE)',
            'icon': 'fas fa-exchange-alt',
            'description': 'Events where one exon is included while the other is skipped, but never both in the same transcript.'
        },
        {
            'id': 'ri',
            'code': 'RI',
            'title': 'Retained Intron (RI)',
            'icon': 'fas fa-pause',
            'description': 'Events where an intron is retained instead of being spliced out.'
        },
        {
            'id': 'mse',
            'code': 'MSE',
            'title': 'Multiple Skipped Exons (MSE)',
            'icon': 'fas fa-layer-group',
            'description': 'Events involving the coordinated skipping of multiple consecutive exons.'
        },
        {
            'id': 'afe',
            'code': 'AFE',
            'title': 'Alternative First Exon (AFE)',
            'icon': 'fas fa-play',
            'description': 'Events involving the use of alternative first exons in transcript isoforms.'
        },
        {
            'id': 'ale',
            'code': 'ALE',
            'title': 'Alternative Last Exon (ALE)',
            'icon': 'fas fa-stop',
            'description': 'Events involving the use of alternative last exons in transcript isoforms.'
        }
    ]
