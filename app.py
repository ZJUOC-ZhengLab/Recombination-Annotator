import pathlib
BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()

import io
import os
import base64

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

import sqlite3
from sqlalchemy.sql import select
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Table, create_engine, Column, String, Integer, MetaData, Index, ForeignKey

from flask_login import LoginManager, login_user, current_user, logout_user, UserMixin
from werkzeug.security import generate_password_hash, check_password_hash

import dash
from dash import dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash import html, dash_table
from dash import no_update


config = dict(
    scrollZoom=True,
    modeBarButtonsToRemove=['lasso2d', 'select2d']
)

conn = sqlite3.connect("annotator.sqlite")
engine = create_engine('sqlite:///annotator.sqlite')
metadata = MetaData()

user_tbl = Table(
    "users",
    metadata,
    Column("id", String(50), primary_key=True),
    Column("username", String(100), unique=True, nullable=False),
    Column("password", String(100), nullable=False)
)

annotation_tbl = Table(
    "annotations",
    metadata,
    Column("id", Integer, primary_key=True, autoincrement=True),
    Column("strain", String(20), nullable=False),
    Column("chrom", Integer, nullable=False),
    Column("event", String(50), nullable=False),
    Column("loh", String(20), nullable=False),
    Column("transition_label", String(20), nullable=False),
    Column("bd_left", Integer, nullable=False),
    Column("bd_right", Integer, nullable=False),
    Column("user_id", String(50), ForeignKey('users.id'), nullable=False)
)

Index("event_group", annotation_tbl.c.strain, annotation_tbl.c.chrom, annotation_tbl.c.event, annotation_tbl.c.loh, annotation_tbl.c.transition_label)

chrom2id = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16}
id2chrom = {v: k for k, v in chrom2id.items()}

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    meta_tags=[{"name": "recombination_annotator", "content": "width=device-width, initial-scale=1"}]
)
server = app.server
app.title = "Recombination Annotator"
app.config.suppress_callback_exceptions = True

server.config.update(
    SECRET_KEY=b'\xc6\x89\xb1\xd5\x03)\x15S\xe0\xdf"\x8a\x11\xf7B\xc5D\x89N\xea\x86\xd2\x8ac',
    SQLALCHEMY_DATABASE_URI='sqlite:///annotator.sqlite',
    SQLALCHEMY_TRACK_MODIFICATIONS=False
)
db = SQLAlchemy()
db.init_app(server)

class Users(db.Model):

    id = db.Column(db.String(50), primary_key=True)
    username = db.Column(db.String(100), unique=True, nullable=False)
    password = db.Column(db.String(100), nullable=False)

    def __init__(self, id, username, password):

        self.id = id
        self.username = username
        self.password = password


login_manager = LoginManager()
login_manager.init_app(server)
login_manager.login_view = '/login'
class Users2(UserMixin, Users):
    pass

@login_manager.user_loader
def load_user(user_id):
    return Users2.query.get(user_id)

login_layout = html.Div([
    # dcc.Location(id="login_url", refresh=False),
    dbc.Container(
        [
            dbc.Row(
                dbc.Col(
                    dbc.Card(
                        [
                            html.Div(
                                dbc.Input(id='login-username',placeholder='Username', style={'height': "34px", "fontSize": "15px"}),
                                style={'margin-top':'10px',}
                            ),
                            html.Div(
                                dbc.Input(id='login-password',placeholder='Password',type='password', style={'height': "34px", "fontSize": "15px"}),
                                style={'margin-top':'10px'}
                            ),
                            html.Div(
                                [
                                    dbc.Button('Log in', id='login-button',color='primary', style={"fontSize": "15px"}),
                                ],
                                className="d-grid gap-2 col-6 mx-auto",
                                style={'margin-top':'10px'}
                            ),
                            html.Br(),
                            html.Div(id='login-alert')
                        ],
                        body=True,
                    ),
                    width=5,
                    style={'margin-top':'30px'}
                ),
                justify='center'
            )
        ]
    )
])


desc_card = html.Div(
    id="description-card",
    children=[
        html.H5("Recombination Annotator"),
        html.H3("Welcome to the Genomic Alteration Annotation Platform"),
        html.Div(
            id="intro",
            children="In the platform, you can visualize sequence data, annotate events, and effectively submit annotation to database.",
        ),
    ],
)


event_type_input = dbc.Row([
    dbc.Label("Event Type", width=4, html_for="event_type_input"),
    dbc.Col(
        dcc.Dropdown(
            id="event_type_input", placeholder="Select event type",
            options={
                "CON": "CON",
                "CO/BIR": "CO/BIR",
                "CON/CO": "CON/CO",
                "COM/CON": "COM/CON",
                "interDEL": "interDEL",
                "terDEL": "terDEL",
                "interDUP": "interDUP",
                "terDUP": "terDUP"
            },
            value="Unknown"
        ),
    ),
], style={"margin-top": "10px"})

LOH_class_input = dbc.Row([
    dbc.Label("LOH Class", width=4, html_for="loh_class_input"),
    dbc.Col(
        dbc.Input(type="text", id="loh_class_input", placeholder="Input LOH class", style={'height': "34px", "fontSize": "15px"}), width=8
    ),
], style={"margin-top": "10px"})

boundary_input = html.Div([
    dbc.Label("Transition Boundary"),
    dbc.Row([
        dbc.Col([
            dbc.Row([
                dbc.Label("Left", width=3, html_for="boundary-left"),
                dbc.Col(
                    dbc.Input(type="text", id="boundary-left", placeholder="0", style={'height': "34px", "fontSize": "15px"}), width=9
                ),
            ])
        ], width=6),
        dbc.Col([
            dbc.Row([
                dbc.Label("Right", width=3, html_for="boundary-right"),
                dbc.Col(
                    dbc.Input(type="text", id="boundary-right", placeholder="0", style={'height': "34px", "fontSize": "15px"}), width=9
                ),
            ])
        ], width=6)
    ])
])

transition_label_input =  dbc.Row([
    dbc.Label("Transition Label", width=5, html_for="transition_label_input"),
    dbc.Col(
        dbc.Input(type="text", id="transition_label_input", placeholder="Input transition label", style={'height': "34px", "fontSize": "15px"}), width=7
    ),
], style={"margin-top": "10px"})


control_card = html.Div(
    id="control-card",
    children=[
        html.P("Upload sequence data"),
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select File')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
            },
        ),
        html.Br(),
        html.P("Selected Point"),
        html.Div(className='row', children=[
            html.Div([
                html.Pre(id='click-data', style = {
                    'border': 'thin lightgrey solid',
                    'overflowX': 'scroll'
                }),
                dbc.FormText("You could set the point as the boundary of transition."),
                dbc.Row([
                    dbc.Col(dbc.Button("Left", id="boundary-left-btn", color="success", style={"fontSize": "15px"}), className="d-grid gap-2", style={'height': '34px'}),
                    dbc.Col(dbc.Button("Right", id="boundary-right-btn", color="danger", style={"fontSize": "15px"}), className="d-grid gap-2", style={'height': '34px'}),
                ]),
            ]),
        ]),

        html.Hr(),
        html.P("Submit Event"),
        boundary_input,
        event_type_input,
        LOH_class_input,
        transition_label_input,
        html.Div(
            children=[
                dbc.Button("Submit", id="submit-btn", color="primary", style={'height': "34px", "fontSize": "15px"}),
                dbc.Alert(
                    "Submit success!",
                    id="success-msg",
                    is_open=False,
                    duration=1000,
                ),
            ],
            className="d-grid gap-2 col-6 mx-auto",
            style={"margin-top": "10px"}
        )
    ]
)


side_panel = html.Div(
    id="left-column",
    className="three columns",
    children=[
        desc_card,
        control_card
    ]
)


loading_spinner = html.Div([
    dbc.Spinner(html.Div(id="loading-output")),
])


search_tab_content = html.Div(
    children=[dbc.Row([
        dbc.Col([
            dbc.Label("Strain", width=4, html_for="search_strain_input"),
            dbc.Input(type="text", id="search_strain_input", placeholder="Input strain",  style={'height': "34px", "font-size": "15px"})
        ], width=3),
        dbc.Col([
            dbc.Label("Chromosome", width=4, html_for="search_chrom_input"),
            dbc.Input(type="text", id="search_chrom_input", placeholder="Input Chromosome",  style={'height': "34px", "font-size": "15px"})
        ], width=3),
        dbc.Col([
            dbc.Label("Event type", width=4, html_for="search_event_input"),
            dcc.Dropdown(
                id="search_event_input", placeholder="Select event type",
                options={
                    "CON": "CON",
                    "CO/BIR": "CO/BIR",
                    "CON/CO": "CON/CO",
                    "COM/CON": "COM/CON",
                    "interDEL": "interDEL",
                    "terDEL": "terDEL",
                    "interDUP": "interDUP",
                    "terDUP": "terDUP"
                },
            ),
        ], width=4),
        dbc.Col([
            html.Div(
                [   
                    dcc.Download(id='data-download'),
                    dbc.Button("Seach", id="search-btn", color="primary", style={"font-size": "15px"}),
                    dbc.Button("Export", id="export-btn", color="secondary", style={"font-size": "15px"}),
                ],
                className="d-grid gap-2",
                style={"margin-top": '5px'}
            )
        ], width=2)
    ]),
    html.Br(),
    html.Div(
        children=[
            dash_table.DataTable(
                data=[],
                columns=[
                    {'name': 'ID', 'id': 'id'},
                    {'name': 'Strain', 'id': 'strain'},
                    {'name': 'Chromosome', 'id': 'chrom'},
                    {'name': 'Event type', 'id': 'event_type'},
                    {'name': 'LOH class', 'id': 'loh'},
                    {'name': 'Transition label', 'id': 'transition_label'},
                    {'name': 'Left', 'id': 'left'},
                    {'name': 'Right', 'id': 'right'}
                ],
                page_size=20,
                style_header={
                    'border': '1px solid #ced4da',
                    'color': '#6c757d',
                    'fontWeight': 'bold',
                    'textAlign': 'center',
                    'backgroundColor': 'white',
                    'font-family': 'Arial'
                },
                style_cell={
                    'border': '1px solid #ced4da',
                    'color': '#6c757d',
                    'textAlign': 'center',
                    'font-family': 'Arial'
                },
                id='search-tbl',
            )
        ]
    ),
])

del_tab_content = html.Div(
    children=[
        dbc.Row([
            dbc.Col([
                dbc.Label("Record Id", width=4, html_for="delete_id_input"),
                dbc.Input(type="text", id="delete_id_input", placeholder="Input record id",  style={'height': "34px", "font-size": "15px"})
            ], width=3),
        ]),
        dbc.Row([
            dbc.Col(
                html.Div(
                    [   
                        dbc.Button("Delete", id="delete-btn", color="danger", style={"font-size": "15px"}),
                        dbc.Alert(
                            "Delete success!",
                            id="delete-msg",
                            is_open=False,
                            duration=1000,
                        ),
                        dbc.Alert(
                            "Try to delete invalid record.",
                            id="delete-fail-msg",
                            color="danger",
                            is_open=False,
                            duration=1000
                        )
                    ],
                    className="d-grid gap-2",
                    style={"margin-top": '5px'}
                ), 
                width=2
            ),
            dbc.Col(
                html.Div(
                    [   
                        dbc.Button("Delete All", id="delete-all-btn", color="dark", style={"font-size": "15px"}),
                        dbc.Alert(
                            "Delete success!",
                            id="delete-all-msg",
                            is_open=False,
                            duration=1000,
                        ),
                    ],
                    className="d-grid gap-2",
                    style={"margin-top": '5px'}
                ), 
                width=2
            )
        ])
    ]
)

advanced_export_tab_content = html.Div(
    children=[
        html.H4("Advanced Export", style={"margin-bottom": "1rem", "margin-top": "1.5rem"}),
        dbc.FormText("Enter strain IDs. Separate IDs by whitespace (space, newline)"),
        dbc.Textarea(
            id="strain_text_area",
            className="mb-3",
            placeholder="WY38#20-1 WY66#30-11\nWY103#15-5",
            style={'height': "100px", "fontSize": "15px"}
        ),
        html.Div(
            children=[
                dcc.Download(id='data-download-2'),
                dbc.Button(
                    "Export IDs",
                    id='ad-export-btn',
                    disabled=True,
                    style={"font-size": "15px", "display": "inline-block"}
                ),
            ],
            style={'text-align': "right"}
        )
    ]
)

tabs = dbc.Tabs(
    [
        dbc.Tab(search_tab_content, label="Search"),
        dbc.Tab(del_tab_content, label="Delete"),
        dbc.Tab(advanced_export_tab_content, label="Export")
    ]
)


fig_panel = html.Div(
    id="right-column",
    className="nine columns",
    children=[   

        html.Div(
            id='data-visual',
            children=[
                html.B("Data Visualization"),
                html.Hr(),
                html.Div([
                    dcc.Store(id='data-memory'),
                    dcc.Store(id='strain-memory'),
                    dcc.Dropdown(
                        id="chrom-input",
                        placeholder="Select chromosome"
                    )
                ]),
                loading_spinner,
                html.Hr(),
                html.Div([
                    dbc.Spinner(html.Div(id="loading-graph"))
                ])
            ]
        ),

        html.Div(
            id='annotation-records',
            children=[
                html.B("Annotation Records"),
                html.Hr(),
                tabs
            ]
        ),
    ]
)


navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src=app.get_asset_url('dna.png'), height="30px", style={"margin-left": "5%"})),
                        dbc.Col(dbc.NavbarBrand("Recombination Annotator", className="ms-2", style={'font-size': "15px"})),
                    ],
                    align="center",
                    className="g-0",
                ),
                href="#",
                style={"textDecoration": "none"},
            ),
            dbc.Row(
                [
                    dbc.Col(dbc.Button("Log out", id='log-out-btn', href='/logout', outline=True, color='light', style={"font-size": "15px",})),
                ],
                className="g-0 flex-nowrap",
                align="center",
            )
        ],
        fluid=True
    ),
    color="dark",
    dark=True,
    style={
        'height': "5.5rem",
        "padding": "0rem 2rem 0rem",
        "align-items": "center"
    }
)


home_layout = html.Div([
    navbar,
    side_panel,
    fig_panel,
])


app.layout = html.Div([
    dcc.Location(id="url", refresh=True),
    html.Div(id="app-container"),
])


def gen_chrom_dropdown(content):

    content_type, content_string = content.split(',')

    decoded = base64.b64decode(content_string)
    try:
        df = pd.read_table(io.StringIO(decoded.decode('utf-8')), sep=' ')
        df.columns = ['chrom', 'pos', 'w303', 'yjm']
        df['chrom'] = df['chrom'].apply(lambda x: x.replace('chr', ''))
        df['pos'] = df['pos'].astype('int64')

    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    chroms = list(df['chrom'].unique())

    chrom_options = {n: n for n in chroms}

    df = df.to_json(orient="split")

    return chrom_options, df


@app.callback(
    Output('app-container', 'children'),
    [Input('url', 'pathname')]
)
def router(pathname):

    if pathname == '/':
        if current_user.is_authenticated:
            return home_layout
        else:
            return login_layout
    elif pathname == '/logout':
        if current_user.is_authenticated:
            logout_user()
            return login_layout
        else:
            return login_layout

    if pathname == '/home':
        if current_user.is_authenticated:
            return home_layout
        else:
            return login_layout
    
    return '404'


@app.callback(
    [Output('url','pathname'),
     Output('login-alert','children')],
    Input('login-button','n_clicks'),
    Input('login-username','value'),
    Input('login-password','value')
)
def login_auth(n_clicks, username, pwd):

    if 'login-button' == dash.ctx.triggered_id:
        user = Users2.query.filter_by(username=username).first()
        if user:
            if check_password_hash(user.password, pwd):
                login_user(user)
                return '/home', dbc.Alert('Sucess!',dismissable=True)

        return no_update, dbc.Alert('Username/Password incorrect.',color='danger',dismissable=True)
    else:
        return no_update, no_update


@app.callback(
    [Output("loading-output", "children"),
    Output("chrom-input", "options"),
    Output("data-memory", "data"),
    Output("strain-memory", "data")],
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    State('upload-data', 'last_modified')
)
def load_output(content, name, date):

    if content is not None:
        chrom_options, data = gen_chrom_dropdown(content)

        name = name[:-5]

        strain_title = html.Div(
            html.H3(f"Strain: {name}", style={"margin-bottom": "1.5rem", "margin-top": "0rem"}),
            style={
                "marginTop": "10px", 
            }
        )
    
        return strain_title, chrom_options, data, name
    else:
        return no_update


@app.callback(
    Output("loading-graph", "children"),
    [Input("chrom-input", "value"),
    Input("data-memory", "data")]
)
def display_graph(chrom, serial_data):

    if chrom is not None:
        chrom_data = pd.read_json(serial_data, orient="split")
        chrom_data = chrom_data[chrom_data['chrom'] == chrom]

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=chrom_data.pos,
                y=chrom_data.w303,
                marker_symbol='circle-open',
                marker=dict(
                    color='rgb(255,0,0)',
                ),
                name='w303'
            )
        )

        fig.add_trace(
            go.Scatter(
                mode='markers',
                x=chrom_data.pos,
                y=chrom_data.yjm,
                marker_symbol='circle-open',
                marker=dict(
                    color='rgb(50,50,255)',
                ),
                name='yjm',
            )
        )

        fig.update_layout(
            height=400,
            showlegend=False,
            margin=dict(l=5, b=5, t=5, r=5),
            plot_bgcolor='white'
        )

        fig.update_yaxes(
            fixedrange=True,
            tick0=-1.0,
            dtick=0.5,
            range=(-0.55, 2.05),
            ticks="outside",
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            zeroline=True,
            zerolinecolor='gray',
            gridcolor='gray', griddash='dot',
            title='Relative coverage',
        )

        fig.update_xaxes(
            ticks="outside",
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            title="Chromosome coordinate",
            exponentformat='none',
            ticksuffix=" Kb",
        )

        return dcc.Graph(
            id="basic-interactions",
            figure=fig,
            config=config,
        )

@app.callback(
    Output('click-data', 'children'),
    Input('basic-interactions', 'clickData'))
def display_click_data(clickData):
    if clickData is None:
        return no_update

    point = clickData['points'][0]
    return f"{point['x']}"


@app.callback(
    Output('boundary-left', 'value'),
    [Input("boundary-left-btn", "n_clicks"),
    Input('click-data', 'children')]
)
def set_left_boundary(n_clicks, selected_point):
    
    if 'boundary-left-btn' == dash.ctx.triggered_id and selected_point is not None:
        return selected_point
    else:
        return no_update
        

@app.callback(
    Output('boundary-right', 'value'),
    [Input("boundary-right-btn", "n_clicks"),
    Input('click-data', 'children')]
)
def set_left_boundary(n_clicks, selected_point):
    
    if 'boundary-right-btn' == dash.ctx.triggered_id and selected_point is not None:
        return selected_point
    else:
        return no_update


@app.callback(
    Output('search-tbl', 'data'),
    Input('search-btn', 'n_clicks'),
    Input('search_strain_input', 'value'),
    Input('search_chrom_input', 'value'),
    Input('search_event_input', "value"),
    State('search-tbl', 'data'),
    State('search-tbl', 'columns'))
def search(n_clicks, strain, chrom, event, rows, columns):

    if chrom is not None:
        if chrom == "":
            chrom = None
        else:
            chrom = chrom2id[chrom]
    
    col_map = {idx: e['id'] for idx, e in enumerate(columns)}
    data = []
    if 'search-btn' == dash.ctx.triggered_id:
        if (strain is None or strain == "") and (chrom is None) and (event is None or event == ""):
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id())
        elif (strain is None or strain == "") and (chrom is None):
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.event == event)
        elif (strain is None or strain == "") and (event is None or event == ""):
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.chrom == chrom)
        elif (event is None or event == "") and (chrom is None):
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.strain == strain)
        elif strain is None or strain == "":
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.event == event, annotation_tbl.c.chrom==chrom)
        elif chrom is None:
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.strain == strain, annotation_tbl.c.event==event)
        elif event is None or event == "":
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.strain == strain, annotation_tbl.c.chrom==chrom)
        else:
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.strain == strain, annotation_tbl.c.chrom==chrom, annotation_tbl.c.event)

        ins = ins.order_by(annotation_tbl.c.strain, annotation_tbl.c.chrom, annotation_tbl.c.event, annotation_tbl.c.loh, annotation_tbl.c.transition_label)
        conn = engine.connect()
        records = conn.execute(ins)
        conn.commit()
        conn.close()

        for rec in records:
            rec_ = {}
            for idx, e in enumerate(rec):
                if idx == 8:
                    continue
                if col_map[idx] == 'chrom':
                    e = id2chrom[e]
                rec_[col_map[idx]] = e
        
            data.append(rec_)

        return data
    else:
        return no_update


@app.callback(
    Output('success-msg', 'is_open'),
    Input('submit-btn', 'n_clicks'),
    Input('strain-memory', 'data'),
    Input('chrom-input', 'value'),
    Input('event_type_input', 'value'),
    Input('loh_class_input', 'value'),
    Input('boundary-left', 'value'),
    Input('boundary-right', 'value'),
    Input('transition_label_input', 'value'),
    State("success-msg", "is_open"),
)
def insert(n_clicks, strain, chrom, event, loh, bd_left, bd_right, tlabel, success_open):

    if 'submit-btn' == dash.ctx.triggered_id:
        
        if strain is not None and chrom is not None and bd_left is not None and bd_right is not None and event is not None and loh is not None and tlabel is not None:
            chrom = chrom2id[chrom]
            
            e = {'strain': strain, 'chrom': chrom, 'event':event,'loh': loh,'transition_label': tlabel, 'bd_left': int(bd_left), 'bd_right': int(bd_right), 'user_id': current_user.get_id()}
            ins = annotation_tbl.insert().values(**e)
            conn = engine.connect()
            conn.execute(ins)
            conn.commit()
            conn.close()

            return True
   
    return False


@app.callback(
    Output('submit-btn', 'disabled'),
    Input('strain-memory', 'data'),
    Input('chrom-input', 'value'),
    Input('event_type_input', 'value'),
    Input('loh_class_input', 'value'),
    Input('boundary-left', 'value'),
    Input('boundary-right', 'value'),
    Input('transition_label_input', 'value'),
)
def form_check(strain, chrom, event, loh, bd_left, bd_right, tlabel):

    if (strain is not None and strain != "") and (chrom is not None and chrom != "") and (
        event is not None and event != "" and event != "Unknown") and (loh is not None and loh != "") and (bd_left is not None or bd_left != "") and (
        bd_right is not None and bd_right != "") and (tlabel is not None and tlabel != ""):
        return False

    else:
        return True


@app.callback(
    Output('data-download', 'data'),
    Input('export-btn', 'n_clicks'),
    State('search-tbl', 'data'),
)
def export(n_clicks, data):
    df = pd.DataFrame.from_dict(data)
    if "export-btn" == dash.ctx.triggered_id:
        return dcc.send_data_frame(df.to_excel, 'annotation.xlsx')

@app.callback(
    [Output('delete-msg', 'is_open'),
    Output('delete-fail-msg', 'is_open')],
    Input('delete_id_input', 'value'),
    Input('delete-btn', 'n_clicks'),
    State("delete-msg", "is_open"),
)
def delete(record_id, n_clicks, is_open):

    if 'delete-btn' == dash.ctx.triggered_id:
        if record_id is not None and record_id != "":
            ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.id == record_id)
            conn = engine.connect()
            recs = list(conn.execute(ins))
            conn.commit()
            conn.close()

            if len(recs) > 0:
                ins = annotation_tbl.delete().where(annotation_tbl.c.id == record_id)
                conn = engine.connect()
                conn.execute(ins)
                conn.commit()
                conn.close()

                return True, False
            else:
                return False, True
        else:
            return False, False
    
    return False, False

@app.callback(
    Output('delete-all-msg', 'is_open'),
    Input('delete-all-btn', 'n_clicks'),
    State("delete-all-msg", "is_open"),
)
def delete_all(n_clicks, is_open):

    if 'delete-all-btn' == dash.ctx.triggered_id:
        
        ins = annotation_tbl.delete().where(annotation_tbl.c.user_id == current_user.get_id())
        conn = engine.connect()
        conn.execute(ins)
        conn.commit()
        conn.close()

        return True
    else:
        return False

@app.callback(
   [Output('ad-export-btn', 'children'),
    Output('ad-export-btn', 'disabled')],
    Input('strain_text_area', 'value')
)
def valid_export_textarea(text):

    if text is None or text == "":
        return "Export IDs", True
    else:
        strains_set = set()
        strain_list = []

        text = text.replace('\n', ' ')
        items = text.split()
        for e in items:
            if e not in strains_set:
                strains_set.add(e)
                strain_list.append(e)

        return f"Export {len(strains_set)} IDs", False

@app.callback(
    Output('data-download-2', 'data'),
    Input('ad-export-btn', 'n_clicks'),
    Input('strain_text_area', 'value')
)
def ad_export(n_clicks, text):

    if 'ad-export-btn' == dash.ctx.triggered_id:
        strains_set = set()
        strain_list = []

        text = text.replace('\n', ' ')
        items = text.split()
        for e in items:
            if e not in strains_set:
                strains_set.add(e)
                strain_list.append(e)

        ins = annotation_tbl.select().where(annotation_tbl.c.user_id == current_user.get_id(), annotation_tbl.c.strain.in_(strain_list))
        ins = ins.order_by(annotation_tbl.c.strain, annotation_tbl.c.chrom, annotation_tbl.c.event, annotation_tbl.c.loh, annotation_tbl.c.transition_label)
        conn = engine.connect()
        records = conn.execute(ins)
        conn.commit()
        conn.close()

        df = pd.DataFrame(records, columns=['ID', 'Strain', 'Chromosome', 'Event type', 'LOH class', 'Transition label', 'Left', 'Right', 'User'])
        df['Chromosome'] = df['Chromosome'].map(id2chrom)
        df = df.drop(columns=['User'])
        return dcc.send_data_frame(df.to_excel, 'annotation.xlsx')
    else:
        return no_update


if __name__ == "__main__":
    app.run_server(debug=False, )
